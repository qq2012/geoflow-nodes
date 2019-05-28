#include "ConnectCandidatePt.h"
#include <queue>
#include <map>


//void ridge::FindTopology(line &smoothLines, const masb::intList &symple_idList, const int_pair_vec &adjacency) {
    /*
    std::cout << "smoothLineid\n";
    for (auto i : symple_idList)
        std::cout << i << " ";
    std::cout << "\n";
    for (auto adj_pair : adjacency) {
        int sheet_id1 = adj_pair.first;
        int sheet_id2 = adj_pair.second;
        if (sheet_id1 >= sheet_id2)
            continue;
        auto it1 = std::find(symple_idList.begin(), symple_idList.end(), sheet_id1);
        auto idx1 = std::distance(symple_idList.begin(), it1);
        auto it2 = std::find(symple_idList.begin(), symple_idList.end(), sheet_id2);
        auto idx2 = std::distance(symple_idList.begin(), it2);

        auto pt1 = smoothLines[idx1][0];
        auto pt2 = smoothLines[idx1].back();

        kdtree2::KDTree* kdtree2 = new kdtree2::KDTree(smoothLines[idx2], true);
        kdtree2->sort_results = true;

        kdtree2::KDTreeResultVector nearest1, nearest2;
        kdtree2->n_nearest(pt1, 1, nearest1);
        kdtree2->n_nearest(pt2, 1, nearest2);


        auto pt3 = smoothLines[idx2][0];
        auto pt4 = smoothLines[idx2].back();

        kdtree2::KDTree* kdtree1 = new kdtree2::KDTree(smoothLines[idx1], true);
        kdtree1->sort_results = true;

        kdtree2::KDTreeResultVector nearest3, nearest4;
        kdtree1->n_nearest(pt3, 1, nearest3);
        kdtree1->n_nearest(pt4, 1, nearest4);

        std::vector<float> dis_vec = { nearest1[0].dis,nearest2[0].dis,nearest3[0].dis,nearest4[0].dis };

        int minElementIndex = std::min_element(dis_vec.begin(), dis_vec.end()) - dis_vec.begin();
        if (minElementIndex == 1) {
            auto i = nearest1[0].idx;
            smoothLines[idx1].insert(smoothLines[idx1].begin(), smoothLines[idx2][i]);
        }
        if (minElementIndex == 2) {
            auto i = nearest2[0].idx;
            smoothLines[idx1].push_back(smoothLines[idx2][i]);
        }
        if (minElementIndex == 3) {
            auto i = nearest3[0].idx;
            smoothLines[idx2].insert(smoothLines[idx2].begin(), smoothLines[idx1][i]);
        }
        if (minElementIndex == 4) {
            auto i = nearest4[0].idx;
            smoothLines[idx2].push_back(smoothLines[idx1][i]);
        }
    }
    */

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
    
}
*/