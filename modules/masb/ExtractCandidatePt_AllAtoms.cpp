#include "ExtractCandidatePt.h"
using namespace masb;
void  masb::allAtoms2Candidates(ExtractCandidatePt_pram & power, MAT &mat, PointList &pointcloud,
    intList &seg_id, PointList &unShrinkingPt, PointList &candidatePt, intList &candidatePt_id) {

    kdtree2::KDTree* pc_kdtree;
    pc_kdtree = new kdtree2::KDTree(pointcloud, true);
    pc_kdtree->sort_results = true;

    kdtree2::KDTree* unShrinkingPt_kdtree;
    unShrinkingPt_kdtree = new kdtree2::KDTree(unShrinkingPt, true);
    unShrinkingPt_kdtree->sort_results = true;

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

        auto cur_sheet_kdtree = new kdtree2::KDTree(cur_sheet_mat.atom, true);
        cur_sheet_kdtree->sort_results = true;
        for (int i = 0; i < cur_sheet_mat.atom.size(); ++i) {
            auto pt = cur_sheet_mat.atom[i];
            kdtree2::KDTreeResultVector neighbours_bis_avg;
            cur_sheet_kdtree->n_nearest(pt, power.bis_avg_knn, neighbours_bis_avg);
            float r_sum = 0;
            for (auto &n : neighbours_bis_avg) {
                r_sum += cur_sheet_mat.radius[n.idx];
            }
            Vector bis_avg = Vector(0.0, 0.0, 0.0);
            for (auto &n : neighbours_bis_avg) {
                bis_avg = bis_avg + (cur_sheet_mat.radius[n.idx] / r_sum) * cur_sheet_mat.bisector[n.idx];
            }
            masb::Point xyz = pt + bis_avg * cur_sheet_mat.radius[i];

            kdtree2::KDTreeResultVector neighbours2PointCloud;
            //pc_kdtree->n_nearest(this->can_pt_r_bisector_avg[i], 1, neighbours);
            pc_kdtree->n_nearest(pt, 1, neighbours2PointCloud);
            auto nearestPt = pointcloud[neighbours2PointCloud[0].idx];

            Point candidate(xyz[0], xyz[1], nearestPt[2]);

            kdtree2::KDTreeResultVector neighbours2unShrinking;
            unShrinkingPt_kdtree->n_nearest(candidate, 1, neighbours2unShrinking);
            if (neighbours2unShrinking[0].dis > power.unshrinkingDist) {
                //if (neighbours2unShrinking[0].dis > 0) {
                candidatePt.push_back(candidate);
                candidatePt_id.push_back(cur_id);
            }
        }

    }
}