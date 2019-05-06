#include "ExtractCandidatePt.h"
using namespace masb;

inline bool ExtractCandidatePt::validateCandidate(float deviationAng_thres, Vector &vec1, Vector &vec2) {
    return vec1 * vec2 > deviationAng_thres;
}

void ExtractCandidatePt::processing(ExtractCandidatePt_pram & power, mat_data &madata,
    ma_Geometry &maGeometry, Sheet_idx_List &sheets, PointCloud &PointCloud) {

    std::cout << "cos(deviationAng_thres) is " << power.deviationAng_thres << std::endl;

    for (int id = 0; id < sheets.size();++id){
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
            while (i_next != -1 && visitFlag[i_next] == 0){
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
            this->seg_id.push_back(id);// seg_id;
            this->direction.push_back(masheetGeometry.ma_normal[idx]);

            //cp only works for interior have bug for exterior MAT
            auto ground_pt_id = a_sheet[idx];
            this->cp.push_back(PointCloud.normals[ground_pt_id]);
            Point q = PointCloud.coords[maSheetData.ma_qidx[idx]];
            auto vcq = q - pt;
            vcq.normalize();
            this->cq.push_back(vcq);

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
            this->can_pt_bisector_avg.push_back(pt + bis_avg * r);

        }
        //this->candidate_cos.push_back(candidatept_cos);
        //this->candidate_r.push_back(candidatept_r);
        //this->candidate_dir.push_back(candidatept_dir);
        //this->candidate_size += candidatept_r.size();
    }

}


