#include "ExtractCandidatePt.h"
using namespace masb;

inline bool ExtractCandidatePt::validateCandidate(float deviationAng_thres, Vector &vec1, Vector &vec2) {
    return vec1 * vec2 > deviationAng_thres;
}

void ExtractCandidatePt::processing_old(ExtractCandidatePt_pram & power, mat_data &madata,
    ma_Geometry &maGeometry, Sheet_idx_List &sheets, PointCloud &PointCloud) {

    std::cout << "cos(deviationAng_thres) is " << power.deviationAng_thres << std::endl;

    for (int id = 0; id < sheets.size(); ++id) {
        //for(auto &a_sheet : sheets) {
        auto a_sheet = sheets[id];
        masb::idx_filter filter;
        mat_data maSheetData;
        ma_Geometry masheetGeometry;
        filter.processing(madata, maGeometry, a_sheet, maSheetData, masheetGeometry);
        //intList ma_pidx; ma_pidx = sheet_idx;
        if (maSheetData.kdtree_ma_coords == NULL)
            maSheetData.kdtree_ma_coords = new kdtree2::KDTree(maSheetData.ma_coords, true);
        maSheetData.kdtree_ma_coords->sort_results = true;

        intList visitFlag; visitFlag.resize(a_sheet.size(), 0);
        size_tList candidatept_idx;
        int i_next = -1;
        for (int i = 0; i < maSheetData.ma_coords.size(); i++) {
            i_next = i;
            while (i_next != -1 && visitFlag[i_next] == 0) {
                visitFlag[i_next] == 1;
                kdtree2::KDTreeResultVector neighbours;
                Point pt = maSheetData.ma_coords[i_next];
                maSheetData.kdtree_ma_coords->r_nearest(pt, power.SearchRadius, neighbours);
                if (neighbours.size() == 1) {
                    i_next = -1;
                }
                else {
                    float r = maSheetData.ma_radius[i_next];
                    std::vector<size_t> smallRindL;
                    for (auto &candidate : neighbours) {
                        if (maSheetData.ma_radius[candidate.idx] < r) {
                            smallRindL.push_back(candidate.idx);
                        }
                    }
                    if (smallRindL.size() == 0) {
                        candidatept_idx.push_back(i_next);
                        i_next = -1;
                    }
                    else {
                        int con = 0;
                        for (auto idx : smallRindL) {
                            Vector vec_pt_bis = masheetGeometry.ma_bisector[i_next];
                            Point candidate_pt = maSheetData.ma_coords[idx];
                            auto vec_ptCand = candidate_pt - pt;
                            vec_ptCand.normalize();
                            if (validateCandidate(power.deviationAng_thres, vec_pt_bis, vec_ptCand)) {
                                i_next = idx;
                                break;
                            }
                            else {
                                con++;
                            }
                        }
                        if (con == smallRindL.size()) {
                            candidatept_idx.push_back(i_next);
                            i_next = -1;
                        }
                    }
                }
            }
        }
        for (auto idx : candidatept_idx) {
            auto pt = maSheetData.ma_coords[idx];
            auto r = maSheetData.ma_radius[idx];
            auto bisec = masheetGeometry.ma_bisector[idx];
            auto SepAng = masheetGeometry.ma_SeperationAng[idx];
            //auto vec = bisec * r;
            //auto can = pt + bisec * r;
            //candidatept_r.push_back(pt + bisec * r);
            //auto can = pt + bisec * (r / cos(SepAng / 2));
            //candidatept_cos.push_back(pt + bisec * (r / cos(SepAng / 2)));
            //candidatept_dir.push_back(masheetGeometry.ma_normal[idx]);
            this->can_pt_r.push_back(pt + bisec * r);
            this->can_pt_cos.push_back(pt + bisec * (r / cos(SepAng / 2)));
            //PointList can_pt_bisector_avg;
            //this->seg_id.push_back(id);// seg_id;
            //this->direction.push_back(masheetGeometry.ma_normal[idx]);

            //cp only works for interior have bug for exterior MAT
            auto ground_pt_id = a_sheet[idx];
            //this->cp.push_back(PointCloud.normals[ground_pt_id]);
            Point q = PointCloud.coords[maSheetData.ma_qidx[idx]];
            auto vcq = q - pt;
            vcq.normalize();
            //this->cq.push_back(vcq);

            kdtree2::KDTreeResultVector neighbours_bis_avg;
            maSheetData.kdtree_ma_coords->n_nearest(pt, 20, neighbours_bis_avg);
            float r_sum = 0;
            for (auto &n : neighbours_bis_avg) {
                r_sum += maSheetData.ma_radius[n.idx];
            }
            Vector bis_avg = Vector(0.0, 0.0, 0.0);
            for (auto &n : neighbours_bis_avg) {
                bis_avg = bis_avg + (maSheetData.ma_radius[n.idx] / r_sum) * masheetGeometry.ma_bisector[n.idx];
            }
            this->can_pt_r_bisector_avg.push_back(pt + bis_avg * r);

        }
        //this->candidate_cos.push_back(candidatept_cos);
        //this->candidate_r.push_back(candidatept_r);
        //this->candidate_dir.push_back(candidatept_dir);
        //this->candidate_size += candidatept_r.size();
    }
}

void ExtractCandidatePt::EdgeBallDetection(ExtractCandidatePt_pram & power, MAT &mat, int cur_id) {
    auto kdtree = new kdtree2::KDTree(mat.atom, true);
    kdtree->sort_results = true;
    for (int idx = 0; idx < mat.matSize; ++idx) {
        kdtree2::KDTreeResultVector neighbours;
        Point pt = mat.atom[idx];
        kdtree->r_nearest(pt, power.SearchRadius, neighbours);
        if (neighbours.size() == 1)
            continue;
        for (int i = neighbours.size() - 1; i >= 0; i--) {
            if (mat.radius[neighbours[0].idx] > mat.radius[idx])
                neighbours.erase(neighbours.begin() + i);
        }
        int counter = 0;
        for (auto &n : neighbours) {
            Vector pt_bisector = mat.bisector[idx];
            Point neighbour = mat.atom[n.idx];
            auto vec_pt2neigh = neighbour - pt;
            vec_pt2neigh.normalize();
            if (validateCandidate(power.deviationAng_thres, pt_bisector, vec_pt2neigh))
                break;
            else
                counter++;
        }
        if (counter == neighbours.size()) {
            this->edgeBalls.atom.push_back(mat.atom[idx]);
            auto r = mat.radius[idx];
            this->edgeBalls.radius.push_back(r);
            //this->edgeBalls.oundPoint;
            this->edgeBalls.fp.push_back(mat.fp[idx]);
            this->edgeBalls.fq.push_back(mat.fq[idx]);
            //this->edgeBalls.p_reverse_norm;
            this->edgeBalls.sp.push_back(mat.sp[idx]);
            this->edgeBalls.sq.push_back(mat.sq[idx]);
            auto bisec = mat.bisector[idx];
            this->edgeBalls.bisector.push_back(bisec);
            auto SepAng = mat.seperationAng[idx];
            this->edgeBalls.seperationAng.push_back(SepAng);
            this->edgeBalls.ma_direction.push_back(mat.ma_direction[idx]);
            this->edgeBall_id.push_back(cur_id);
            ///////////////////////////////////////////////////////////////
            //                     EdgeBall2CandidatePt
            ///////////////////////////////////////////////////////////////
            this->can_pt_r.push_back(pt + bisec * r);
            this->can_pt_cos.push_back(pt + bisec * (r / cos(SepAng / 2)));

            kdtree2::KDTreeResultVector neighbours_bis_avg;
            kdtree->n_nearest(pt, power.bis_avg_knn, neighbours_bis_avg);
            float r_sum = 0;
            for (auto &n : neighbours_bis_avg) {
                r_sum += mat.radius[n.idx];
            }
            Vector bis_avg = Vector(0.0, 0.0, 0.0);
            for (auto &n : neighbours_bis_avg) {
                bis_avg = bis_avg + (mat.radius[n.idx] / r_sum) * mat.bisector[n.idx];
            }
            this->can_pt_r_bisector_avg.push_back(pt + bis_avg * r);
        }
    }
}


void ExtractCandidatePt::filterCandidatePt(ExtractCandidatePt_pram & power, PointList &pointcloud) {
    kdtree2::KDTree* pc_kdtree;
    pc_kdtree = new kdtree2::KDTree(pointcloud, true);
    pc_kdtree->sort_results = true;
    std::vector<float> dis_list; dis_list.reserve(edgeBalls.matSize);
    for (int i = 0; i < edgeBalls.matSize; ++i) {
        kdtree2::KDTreeResultVector neighbours;
        pc_kdtree->n_nearest(this->can_pt_r_bisector_avg[i], 1, neighbours);
        dis_list.push_back(neighbours[0].dis);
        if (neighbours[0].dis > power.filterDistance) {
            filter.push_back(0);
        }
        else {
            filter.push_back(1);
            //this->remainingPt.push_back(this->can_pt_r_bisector_avg[i]);
            this->candidate_id.push_back(edgeBall_id[i]);
            ///////////////////////////////////////////////////////////////
            //                     NN spatialInterpolation
            ///////////////////////////////////////////////////////////////
            auto xyz = this->can_pt_r_bisector_avg[i];
            auto z_NN = pointcloud[neighbours[0].idx][2];
            this->candidate_pt.push_back({ xyz[0],xyz[1], z_NN });
        }
    }
}
/*
void ExtractCandidatePt::spatialInterpolation(PointList &pointcloud) {

}
*/
void ExtractCandidatePt::processing(ExtractCandidatePt_pram & power, MAT &mat, PointList &pointcloud, intList &seg_id) {
    std::cout << "cos(deviationAng_thres) is " << power.deviationAng_thres << std::endl;
    //auto it = max_element(std::begin(seg_id), std::end(seg_id));
    auto max_id = *max_element(seg_id.begin(), seg_id.end());
    for (int cur_id = 1; cur_id <= max_id; ++cur_id) {
        MAT cur_sheet_mat;
        cur_sheet_mat.isInterior = mat.isInterior;
        for (int idx = 0; idx < seg_id.size(); ++idx) {
            if (seg_id[idx] == cur_id) {
                cur_sheet_mat.atom.push_back(mat.atom[idx]);
                cur_sheet_mat.radius.push_back(mat.radius[idx]);
                //unshrikingGroundPoint;
                cur_sheet_mat.fp.push_back(mat.fp[idx]);
                cur_sheet_mat.fq.push_back(mat.fq[idx]);
                //VectorList sp_reverse_norm;
                cur_sheet_mat.sp.push_back(mat.sp[idx]);
                cur_sheet_mat.sq.push_back(mat.sq[idx]);
                cur_sheet_mat.bisector.push_back(mat.bisector[idx]);
                cur_sheet_mat.seperationAng.push_back(mat.seperationAng[idx]);
                cur_sheet_mat.ma_direction.push_back(mat.ma_direction[idx]);
            }
        }
        cur_sheet_mat.matSize = cur_sheet_mat.atom.size();
        EdgeBallDetection(power, cur_sheet_mat, cur_id);
    }
    this->edgeBalls.matSize = edgeBalls.atom.size();
    filterCandidatePt(power, pointcloud);
    //spatialInterpolation(pointcloud);
}