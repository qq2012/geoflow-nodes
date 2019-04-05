#include "ConnectCandidatePt.h"

#include <unsupported\Eigen\Splines>
//#include "C:\dev\Eigen3\include\eigen3\unsupported\Eigen\Splines"
//#include "C:\dev\Eigen3\include\eigen3\Eigen\Core"

#include <iostream>
#include <fstream>
#include <map>

using namespace Eigen;

void ridge::connectCandidatePt8Spline(masb::PointList &PointCloud, masb::PointList &candidate,
    masb::intList &seg_id, masb::intList &filter, ridge::line &line_segmentList, masb::intList &idList) {
    float filter_thresh = 5.0;

    kdtree2::KDTree* pc_kdtree;
    pc_kdtree = new kdtree2::KDTree(PointCloud, true);
    pc_kdtree->sort_results = true;
    std::vector<float> dis_list; dis_list.reserve(candidate.size());
    for (int i = candidate.size() - 1; i >= 0; i--) {
        kdtree2::KDTreeResultVector neighbours;
        pc_kdtree->n_nearest(candidate[i], 1, neighbours);
        dis_list.push_back(neighbours[0].dis);
        if (neighbours[0].dis > filter_thresh) {
            filter.push_back(0);
            candidate.erase(candidate.begin() + i);
            seg_id.erase(seg_id.begin() + i);
        }
        else
            filter.push_back(1);
    }
    //dis_list;

    std::map<int, int> seg_frequency;
    for (int i : seg_id)
        ++seg_frequency[i];
    //for (const auto& e : seg_frequency)
    //    std::cout << "Element " << e.first 
    //     << " encountered " << e.second << " times\n";

    for (const auto& e : seg_frequency) {
        auto cur_sheet = e.first;
        auto cur_size = e.second;
        masb::PointList cur_candidate;
        cur_candidate.reserve(cur_size);
        for (int i = 0; i < seg_id.size(); i++) {
            if (seg_id[i] == cur_sheet) {
                cur_candidate.push_back(candidate[i]);
            }
        }
        if (cur_candidate.size() != cur_size)
            std::cout << "Error" << std::endl;

        typedef Spline<float, 3> Spline3d;
        typedef Spline3d::PointType PointType;
        typedef Spline3d::KnotVectorType KnotVectorType;
        typedef Spline3d::ControlPointVectorType ControlPointVectorType;

        Eigen::MatrixXf matrix_pt(3, cur_size);
        int col = 0;
        for (auto& p : cur_candidate) {
            matrix_pt.col(col++) << p[0], p[1], p[2];
            //col++;
        }
        //std::cout << matrix_pt;
        ControlPointVectorType ctrls_points;// = ControlPointVectorType::Random(3, 100);
        ctrls_points = matrix_pt;
        
        const Spline3d spline = SplineFitting<Spline3d>::Interpolate(ctrls_points, 3);

        KnotVectorType chord_lengths; // knot parameters
        Eigen::ChordLengths(ctrls_points, chord_lengths);

        for (Eigen::DenseIndex i = 0; i < cur_size; ++i)
        {
            PointType pt = spline(chord_lengths(i));
            PointType ref = matrix_pt.col(i);
            //VERIFY((pt - ref).matrix().norm() < 1e-14);
        }

        auto a_line = Spline3d(chord_lengths, ctrls_points);

    }
}