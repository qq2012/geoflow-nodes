#include "Ma_utility.h"
#include <iostream>

using namespace masb;

intList masb::conditionFilter(masb::MAT &mat, condition_pram & power, PointList &pointcloud,
    PointList &unShrinkingPt, intList &seg_id){

    intList filter;

    kdtree2::KDTree* pc_kdtree;
    pc_kdtree = new kdtree2::KDTree(pointcloud, true);
    pc_kdtree->sort_results = true;

    kdtree2::KDTree* unShrinkingPt_kdtree;
    unShrinkingPt_kdtree = new kdtree2::KDTree(unShrinkingPt, true);
    unShrinkingPt_kdtree->sort_results = true;

    for (int i = 0; i < mat.matSize; ++i) {
        if (seg_id[i] == 0) {
            filter.push_back(0);
            continue;
        }

        auto p = mat.atom[i];
        kdtree2::KDTreeResultVector neighbours;
        pc_kdtree->n_nearest(p, 1, neighbours);

        auto b = mat.bisector[i];
        auto r = mat.radius[i];
        auto c = p + r * b;
        kdtree2::KDTreeResultVector neighbours2;
        unShrinkingPt_kdtree->n_nearest(c, 1, neighbours2);
        //////////////////////////////////////////////////////////////////////////////////////////////
        //     distance2ground, filterEdgeBall8radius and distance2unShrinkingPt
        //////////////////////////////////////////////////////////////////////////////////////////////
        if (neighbours[0].dis < power.filterDistance
            && mat.radius[i] < power.MaxEdgeBallRadius
            && mat.radius[i] > power.MinEdgeBallRadius
            && neighbours2[0].dis > power.unshrinkingDist) {
            filter.push_back(1);
        }
        else {
            filter.push_back(0);
        }
    }
    return filter;
}