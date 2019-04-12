#include "ConnectCandidatePt.h"
#include <queue>


void ridge::FindTopology(line &smoothLines) {
    ////////////////////////////////////////
    float connect_thresh = 40;
    ////////////////////////////////////////
    for (int i = 0; i < smoothLines.size();++i) {
        auto pt1 = smoothLines[i][0];
        auto pt2 = smoothLines[i].back();

        masb::PointList kd_pt;
        for (int j = 0; j < smoothLines.size(); ++j) {
            if (j == i)
                continue;
            //auto a = smoothLines[j];
            kd_pt.insert(kd_pt.end(), smoothLines[j].begin(), smoothLines[j].end());
        }
        kdtree2::KDTree* kdtree = new kdtree2::KDTree(kd_pt, true);
        kdtree->sort_results = true;
        
        kdtree2::KDTreeResultVector nearest1, nearest2;
        kdtree->n_nearest(pt1, 1, nearest1);
        kdtree->n_nearest(pt2, 1, nearest2);

        if (pt1[2]> pt2[2]){
            if (nearest1[0].dis < connect_thresh) {
                auto idx = nearest1[0].idx;
                smoothLines[i].insert(smoothLines[i].begin(), kd_pt[idx]);
            }
        }
        else {
            if (nearest2[0].dis < connect_thresh) {
                auto idx = nearest2[0].idx;
                smoothLines[i].push_back(kd_pt[idx]);
            }
        }
    }
}