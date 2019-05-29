#include "ConnectCandidatePt.h"
#include <queue>
#include <map>
#include <algorithm>

void ridge::FindTopology(PolyineList &polylines, const masb::intList &polylinesID, 
    const int_pair_vec &adjacency){
    
    std::cout << "ridge::FindTopology()\n"
        <<"polylinesID\n";
    for (auto i : polylinesID)
        std::cout << i << " ";
    std::cout << "\n";

    for (auto adj_pair : adjacency) {
        int sheet_id1 = adj_pair.first;//the smaller id
        int sheet_id2 = adj_pair.second;//the bigger id

        if ((std::find(polylinesID.begin(), polylinesID.end(), sheet_id1) == polylinesID.end())
            || (std::find(polylinesID.begin(), polylinesID.end(), sheet_id2) == polylinesID.end()))
            /* polylinesID does not contains sheet_id1 or sheet_id2 */
            continue;

        auto it1 = std::find(polylinesID.begin(), polylinesID.end(), sheet_id1);
        auto idx1 = std::distance(polylinesID.begin(), it1);
        auto it2 = std::find(polylinesID.begin(), polylinesID.end(), sheet_id2);
        auto idx2 = std::distance(polylinesID.begin(), it2);

        auto pt11 = polylines[idx1][0];
        auto pt12 = polylines[idx1].back();

        kdtree2::KDTree* kdtree2 = new kdtree2::KDTree(polylines[idx2], true);
        kdtree2->sort_results = true;

        kdtree2::KDTreeResultVector nearest11, nearest12;
        kdtree2->n_nearest(pt11, 1, nearest11);
        kdtree2->n_nearest(pt12, 1, nearest12);

        auto pt21 = polylines[idx2][0];
        auto pt22 = polylines[idx2].back();

        kdtree2::KDTree* kdtree1 = new kdtree2::KDTree(polylines[idx1], true);
        kdtree1->sort_results = true;

        kdtree2::KDTreeResultVector nearest21, nearest22;
        kdtree1->n_nearest(pt21, 1, nearest21);
        kdtree1->n_nearest(pt22, 1, nearest22);

        std::vector<float> dis_vec = { nearest11[0].dis,nearest12[0].dis,nearest21[0].dis,nearest22[0].dis };

        int minElementIndex = std::min_element(dis_vec.begin(), dis_vec.end()) - dis_vec.begin();
        if (minElementIndex == 0) {
            auto i = nearest11[0].idx;//idx of the nearest point in polyline[idx2]
            polylines[idx1].insert(polylines[idx1].begin(), polylines[idx2][i]);
        }
        if (minElementIndex == 1) {
            auto i = nearest12[0].idx;//idx of the nearest point in polyline[idx2]
            polylines[idx1].push_back(polylines[idx2][i]);
        }
        if (minElementIndex == 2) {
            auto i = nearest21[0].idx;//idx of the nearest point in polyline[idx1]
            polylines[idx2].insert(polylines[idx2].begin(), polylines[idx1][i]);
        }
        if (minElementIndex == 3) {
            auto i = nearest22[0].idx;//idx of the nearest point in polyline[idx1]
            polylines[idx2].push_back(polylines[idx1][i]);
        }
    }
    
    ////////////////////////////////////////
    //float connect_thresh = 80;
    ////////////////////////////////////////
    /*
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
    */
}
