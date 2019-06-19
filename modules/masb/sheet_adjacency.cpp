#include "sheet_adjacency.h"
#include<set>

using namespace masb;
using namespace std;
void ridge::adjacencyProcessing(ridge::adjacency_parameters& power, masb::MAT &mat, masb::PointList &pointCloud,
    masb::Sheet_idx_List &sheets, int_pair_vec &adjacency, masb::PointList &junction, masb::intList &junctionID) {

    if (sheets.size() == 1) {
        std::cout << "ridge::adjacencyProcessing() only one sheet" << std::endl;
        return;
    }
    /*
    PointList pointcloud2d;
    for (auto &pt : pointCloud)
        pointcloud2d.push_back(Point(pt[0], pt[1], 0));
    kdtree2::KDTree* pc2d_kdtree;
    pc2d_kdtree = new kdtree2::KDTree(pointcloud2d, true);
    pc2d_kdtree->sort_results = true;
    */
    kdtree2::KDTree* pc_kdtree;
    pc_kdtree = new kdtree2::KDTree(pointCloud, true);
    pc_kdtree->sort_results = true;

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


            masb::Point JunctionPt4maxZ_ij = masb::Point(0, 0, 0);

            //find points belong to i but have neighbours in j
            std::set<int> ptIdxSet_j;
            for (auto &pt : pt_in_sheet_i) {
                kdtree2::KDTreeResultVector neighbours;
                sheet_j_kdtree->r_nearest(pt, power.searchRadius, neighbours);
                for (auto &n: neighbours) {
                    ptIdxSet_j.insert(n.idx);
                }
            }
            std::set<int>::iterator it_j;
            for (it_j = ptIdxSet_j.begin(); it_j != ptIdxSet_j.end(); it_j++) {
                //auto a = *it_j;
                junction.push_back(pt_in_sheet_j[*it_j]);
                if (power.OnlySurfaceAdjacent &&
                    JunctionPt4maxZ_ij[2] < pt_in_sheet_j[*it_j][2]) {
                    JunctionPt4maxZ_ij = pt_in_sheet_j[*it_j];
                }
            }
            junctionID.insert(junctionID.begin(), ptIdxSet_j.size(), sheet_j_id);

            //find points belong to j but have neighbour in i
            std::set<int> ptIdxSet_i;
            for (auto &pt : pt_in_sheet_j) {
                kdtree2::KDTreeResultVector neighbours;
                sheet_i_kdtree->r_nearest(pt, power.searchRadius, neighbours);
                for (auto &n : neighbours) {
                    ptIdxSet_i.insert(n.idx);
                }
            }
            std::set<int>::iterator it_i;
            for (it_i = ptIdxSet_i.begin(); it_i != ptIdxSet_i.end(); it_i++) {
                //auto a = *it_j;
                junction.push_back(pt_in_sheet_i[*it_i]);
                if (power.OnlySurfaceAdjacent &&
                    JunctionPt4maxZ_ij[2] < pt_in_sheet_i[*it_i][2]) {
                    JunctionPt4maxZ_ij = pt_in_sheet_i[*it_i];
                }
            }
            junctionID.insert(junctionID.begin(), ptIdxSet_i.size(), sheet_i_id);

            if (ptIdxSet_i.size() > power.adjacency_thresh && ptIdxSet_j.size() > power.adjacency_thresh) {
                if (power.OnlySurfaceAdjacent) {
                    kdtree2::KDTreeResultVector neighbours;
                    pc_kdtree->n_nearest(JunctionPt4maxZ_ij, 1, neighbours);
                    std::cout << "the point with max z <---> point cloud distance = " << neighbours[0].dis << std::endl;
                    if (neighbours[0].dis < power.DeepSurf_thresh) {
                        adjacency.push_back(std::make_pair(sheet_i_id, sheet_j_id));
                        std::cout << "sheet " << sheet_i_id << " and " << sheet_j_id << " are surface adjacent"
                            << "sheet " << sheet_i_id << " has " << ptIdxSet_i.size() << " neighbours in " << sheet_j_id
                            << "sheet " << sheet_j_id << " has " << ptIdxSet_j.size() << " neighbours in " << sheet_i_id << std::endl;
                    }
                }
                else {
                    adjacency.push_back(std::make_pair(sheet_i_id, sheet_j_id));
                    std::cout << "sheet " << sheet_i_id << " and " << sheet_j_id << " are adjacent"
                        << "sheet " << sheet_i_id << " has " << ptIdxSet_i.size() << " neighbours in " << sheet_j_id
                        << "sheet " << sheet_j_id << " has " << ptIdxSet_j.size() << " neighbours in " << sheet_i_id << std::endl;
                }
            }            
        }
    }
    if (power.OnlySurfaceAdjacent)
        std::cout << "There are " << adjacency.size() << " pairs of surface adjacent sheets" << std::endl;
    else
        std::cout << "There are " << adjacency.size() << " pairs of adjacent sheets" << std::endl;
}