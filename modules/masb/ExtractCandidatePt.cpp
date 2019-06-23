#include "ExtractCandidatePt.h"

using namespace masb;

inline bool ExtractCandidatePt::validateCandidate(float deviationAng_thres, Vector &vec1, Vector &vec2) {
    return vec1 * vec2 > deviationAng_thres;
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

            //////////////////////////////////////////////////////////////////////////////////////
            //                     EdgeBall2 imtermediate CandidatePt
            //////////////////////////////////////////////////////////////////////////////////////
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


void ExtractCandidatePt::filterCandidatePt(ExtractCandidatePt_pram & power, 
    PointList &pointcloud, PointList &unShrinkingPt) {

    kdtree2::KDTree* pc_kdtree;
    pc_kdtree = new kdtree2::KDTree(pointcloud, true);
    pc_kdtree->sort_results = true;
    std::vector<float> BallTerrainDis_list; 
    BallTerrainDis_list.reserve(edgeBalls.matSize);

    PointList pointcloud2d;
    for (auto &pt : pointcloud)
        pointcloud2d.push_back(Point(pt[0], pt[1], 0));
    kdtree2::KDTree* pc2d_kdtree;
    pc2d_kdtree = new kdtree2::KDTree(pointcloud2d, true);
    pc2d_kdtree->sort_results = true;

    kdtree2::KDTree* unShrinkingPt_kdtree;
    unShrinkingPt_kdtree = new kdtree2::KDTree(unShrinkingPt, true);
    unShrinkingPt_kdtree->sort_results = true;
    std::vector<float> cpt2unShrDis_list;//distance from candidate pt (x,y,z_interpolation) to unshrinking point 
    cpt2unShrDis_list.reserve(edgeBalls.matSize);

    for (int i = 0; i < edgeBalls.matSize; ++i) {

        kdtree2::KDTreeResultVector neighbours;
        //pc_kdtree->n_nearest(this->can_pt_r_bisector_avg[i], 1, neighbours);
        pc_kdtree->n_nearest(this->edgeBalls.atom[i], 1, neighbours);
        BallTerrainDis_list.push_back(neighbours[0].dis);

        ///////////////////////////////////////////////////////////////
        //                     NN spatialInterpolation
        ///////////////////////////////////////////////////////////////
        auto xyz = this->can_pt_r_bisector_avg[i];
        kdtree2::KDTreeResultVector neighbours1;
        pc2d_kdtree->n_nearest(masb::Point(xyz[0], xyz[1], 0), 1, neighbours1);
        auto z_NN = pointcloud[neighbours1[0].idx][2];
        Point candidate(xyz[0], xyz[1], z_NN);

        kdtree2::KDTreeResultVector neighbours2;
        unShrinkingPt_kdtree->n_nearest(candidate, 1, neighbours2);
        //////////////////////////////////////////////////////////////////////////////////////////////
        //     distance2ground, filterEdgeBall8radius and distance2unShrinkingPt
        //////////////////////////////////////////////////////////////////////////////////////////////
        if (neighbours[0].dis < power.filterDistance 
            && edgeBalls.radius[i] < power.MaxEdgeBallRadius
            && edgeBalls.radius[i] > power.MinEdgeBallRadius
            && neighbours2[0].dis > power.unshrinkingDist) {
            filter.push_back(1);
            //this->remainingPt.push_back(this->can_pt_r_bisector_avg[i]);
            this->candidate_id.push_back(edgeBall_id[i]);
            this->candidate_pt.push_back(candidate);
            this->candidate_radius.push_back(this->edgeBalls.radius[i]);
            this->candidate_direction.push_back(this->edgeBalls.ma_direction[i]);
        }
        else {
            filter.push_back(0);
        }
    }
}
/*
void ExtractCandidatePt::spatialInterpolation(PointList &pointcloud) {

}
*/
void ExtractCandidatePt::processing(ExtractCandidatePt_pram & power, MAT &mat, PointList &pointcloud, 
    intList &seg_id, PointList &unShrinkingPt) {
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
    filterCandidatePt(power, pointcloud, unShrinkingPt);
    //spatialInterpolation(pointcloud);
}