#include "ConnectCandidatePt.h"
#include <map>

void ridge::ConnectCandidatePt8PolynomialFitting(masb::PointList &pts, masb::intList &pt_id, masb::PointList &pointCloud,
    float error_thresh, ridge::PolyineList &polylines, const int &order) {
    std::cout << "ridge::ConnectCandidatePt8PolynomialFitting()"
        << "\nerror_thresh = " << error_thresh << std::endl;

    masb::PointList pointcloud2d;
    for (auto &pt : pointCloud)
        pointcloud2d.push_back(masb::Point(pt[0], pt[1], 0));
    kdtree2::KDTree* pc2d_kdtree;
    pc2d_kdtree = new kdtree2::KDTree(pointcloud2d, true);
    pc2d_kdtree->sort_results = true;

    std::map<int, int> seg_frequency;
    for (auto i : pt_id) {
        ++seg_frequency[i];
    }
    /*
    for (const auto& e : seg_frequency)
        std::cout << "Element " << e.first
         << " encountered " << e.second << " times\n";
    */
    for (const auto& e : seg_frequency) {
        auto cur_sheet = e.first;
        auto cur_size = e.second;
        if (cur_sheet == 0 || cur_size < 3)
            continue;

        std::cout << "cur_sheet ID = " << cur_sheet
            << "cur_size " << e.second << std::endl;

        masb::PointList cur_pt;
        cur_pt.reserve(cur_size);
        for (int i = 0; i < pt_id.size(); i++) {
            if (pt_id[i] == cur_sheet) {
                cur_pt.push_back(pts[i]);
            }
        }
        if (cur_pt.size() != cur_size)
            std::cout << "Error" << std::endl;

        // coordinate transform

        masb::PointList fittinyPolyline_xy0;
        float min_error;
        masb::PolynomialFitting(cur_pt, order,fittinyPolyline_xy0, min_error);


        std::cout << " average_min_error = min_error / cur_size = " << min_error / cur_size << std::endl;

        if (min_error / cur_size < error_thresh) {
            masb::PointList fittinyPolyline;
            for (auto &xy0 : fittinyPolyline_xy0) {
                kdtree2::KDTreeResultVector neighbours;
                pc2d_kdtree->n_nearest(xy0, 1, neighbours);
                auto z_NN = pointCloud[neighbours[0].idx][2];
                masb::Point p(xy0[0], xy0[1], z_NN);

                // coordinate transform 
                fittinyPolyline.push_back(p);
            }
            polylines.push_back(fittinyPolyline);
        }
        else
            std::cout << "sheet " << cur_sheet << "does not include" << std::endl;
    }
}