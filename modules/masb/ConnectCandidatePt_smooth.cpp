#include "ConnectCandidatePt.h"

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Point_with_normal_3.h>
#include <CGAL/property_map.h>
#include <CGAL/Shape_detection_3.h>
#include <CGAL/structure_point_set.h>
#include <CGAL/Point_set_3.h>
#include <CGAL/bilateral_smooth_point_set.h>
#include <iostream>
#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel  Kernel;
typedef Kernel::Point_3                                      Point;
typedef Kernel::Vector_3                                     Vector;
typedef std::pair<Point, Vector>                              Point_with_normal;
typedef std::vector<Point_with_normal>                       Pwn_vector;
typedef CGAL::First_of_pair_property_map<Point_with_normal>  Point_map;
typedef CGAL::Second_of_pair_property_map<Point_with_normal> Normal_map;

void ridge::connectCandidatePtSmooth(line &symple_segmentList, line &smoothLine) {
    for (auto &a_line:symple_segmentList) {
        Pwn_vector pwn_vec;
        std::cout << "the size of input path: " << a_line.size() << std::endl;
        Point_with_normal pwn;
        Point point(a_line[0][0], a_line[0][1], a_line[0][2]);
        Vector normal = Vector(a_line[0][0] - a_line[1][0], //nx
                                a_line[0][1] - a_line[1][1], //ny
                                a_line[0][2] - a_line[1][2]);//nz
        pwn_vec.push_back(std::make_pair(point, normal));

        for (int i = 1; i < a_line.size(); ++i) {
            Point_with_normal pwn;
            Point point(a_line[i][0], a_line[i][1], a_line[i][2]);
            Vector normal = Vector(a_line[i-1][0]- a_line[i][0], //nx
                                   a_line[i-1][1]- a_line[i][1], //ny
                                   a_line[i-1][2]- a_line[i][2]);//nz
            pwn_vec.push_back(std::make_pair(point, normal));
        }


        // Algorithm parameters
        int k = 20;                 // size of neighborhood. The bigger the smoother the result will be.
                                     // This value should bigger than 1.
        double sharpness_angle = 45; // control sharpness of the result.
                                     // The bigger the smoother the result will be
                                     // must smaller than 90
        int iter_number = 3;         // number of times the projection is applied

        for (int i = 0; i < iter_number; ++i) {
            CGAL::bilateral_smooth_point_set<CGAL::Sequential_tag>(
                pwn_vec, k,
                CGAL::parameters::point_map(CGAL::First_of_pair_property_map<Point_with_normal>()).
                normal_map(CGAL::Second_of_pair_property_map<Point_with_normal>()).sharpness_angle(sharpness_angle));
        }

        masb::PointList a_smooth_line;
        for (auto &pwn_tmp : pwn_vec) {
            auto pt_cgal = pwn_tmp.first;
            //std::cout << pt[0] << " " << pt[1] << " " << pt[2] << std::endl;
            masb::Point pt = masb::Point(pt_cgal[0], pt_cgal[1], pt_cgal[2]);
            a_smooth_line.push_back(pt);
        }
        //std::cout << std::endl;
        smoothLine.push_back(a_smooth_line);
        std::cout << "the size of smooth path: " << a_smooth_line.size() << std::endl;
    }
}