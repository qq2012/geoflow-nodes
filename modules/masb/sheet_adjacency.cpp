#include "sheet_adjacency.h"
#include<set>

using namespace masb;
using namespace std;
void ridge::adjacencyProcessing(ridge::adjacency_parameters& power, masb::MAT &mat, 
    masb::Sheet_idx_List &sheets, int_pair_vec &adjacency, masb::PointList &junction) {

    if (sheets.size() == 1) {
        std::cout << "ridge::adjacencyProcessing only one sheet" << std::endl;
        return;
    }

    for (int i = 0; i < sheets.size(); ++i) {
        int sheet_i_id = i + 1;
        masb::PointList pt_in_sheet_i;
        for (auto &idx : sheets[i]) {
            pt_in_sheet_i.push_back(mat.atom[idx]);
        }
        kdtree2::KDTree* sheet_i_kdtree;
        sheet_i_kdtree = new kdtree2::KDTree(pt_in_sheet_i, true);
        sheet_i_kdtree->sort_results = true;
        for (int j = i + 1; j < sheets.size(); ++j) {
            int sheet_j_id = j + 1;
            masb::PointList pt_in_sheet_j;
            for (auto &idx : sheets[j]) {
                pt_in_sheet_j.push_back(mat.atom[idx]);
            }
            kdtree2::KDTree *sheet_j_kdtree;
            sheet_j_kdtree = new kdtree2::KDTree(pt_in_sheet_j, true);
            sheet_j_kdtree->sort_results = true;

            //find points belong to i but have neighbours in j
            
            std::set<int> ptIdxSet_j;
            int count_ij = 0;
            for (auto &pt : pt_in_sheet_i) {
                kdtree2::KDTreeResultVector neighbours;
                sheet_j_kdtree->r_nearest(pt, power.searchRadius, neighbours);
                count_ij += neighbours.size();
                for (auto &n: neighbours) {
                    ptIdxSet_j.insert(n.idx);
                }
            }
            std::set<int>::iterator it_j;
            for (it_j = ptIdxSet_j.begin(); it_j != ptIdxSet_j.end(); it_j++) {
                //auto a = *it_j;
                junction.push_back(pt_in_sheet_j[*it_j]);
            }


            //find points belong to j but have neighbour in i
            std::set<int> ptIdxSet_i;
            int count_ji = 0;
            for (auto &pt : pt_in_sheet_j) {
                kdtree2::KDTreeResultVector neighbours;
                sheet_i_kdtree->r_nearest(pt, power.searchRadius, neighbours);
                count_ji += neighbours.size();
                for (auto &n : neighbours) {
                    ptIdxSet_i.insert(n.idx);
                }
            }
            std::set<int>::iterator it_i;
            for (it_i = ptIdxSet_i.begin(); it_i != ptIdxSet_i.end(); it_i++) {
                //auto a = *it_j;
                junction.push_back(pt_in_sheet_i[*it_i]);
            }
            std::cout << "sheet " << sheet_i_id << " has " << count_ij << " neighbours in " << sheet_j_id << std::endl;
            std::cout << "sheet " << sheet_j_id << " has " << count_ji << " neighbours in " << sheet_i_id << std::endl;              

        }
    }
}