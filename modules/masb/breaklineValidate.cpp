#include "breaklineValidate.h"

#include <CGAL/Simple_cartesian.h>
typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_3 Point_3;
typedef Kernel::Segment_3 Segment_3;

void ridge::breaklineValidateProcess(masb::PointList &TestVertices, masb::PointList &ReferenceVertices) {
    /*
    Point_3 p(0, 0, 0);
    Point_3 n1(1, 0, 0), n2(0, 1, 0);
    Segment_3 s(n1, n2);
    auto distance_squared = CGAL::squared_distance(s, p);
    */
    float accu_distance = 0;
    auto ReferenceVertices_kdtree = new kdtree2::KDTree(ReferenceVertices, true);
    ReferenceVertices_kdtree ->sort_results = true;

    int n = ReferenceVertices.size();
    for (auto &v :TestVertices) {
        Point_3 v_cgal(v[0], v[1], v[2]);
        float d = 0;

        kdtree2::KDTreeResultVector neighbours;
        ReferenceVertices_kdtree->n_nearest(v, 1, neighbours);
        if (neighbours[0].idx == 0) {
            auto n1 = ReferenceVertices[0];
            auto n2 = ReferenceVertices[1];
            Point_3 n1_cgal(n1[0],n1[1], n1[2]), n2_cgal(n2[0],n2[1],n2[2]);
            Segment_3 s(n1_cgal, n2_cgal);
            auto distance_squared = CGAL::squared_distance(s, v_cgal);
            d = std::sqrt(distance_squared);
        }
        else if (neighbours[0].idx == n - 1) {
            auto n1 = ReferenceVertices[n - 2];
            auto n2 = ReferenceVertices[n - 1];
            Point_3 n1_cgal(n1[0], n1[1], n1[2]), n2_cgal(n2[0], n2[1], n2[2]);
            Segment_3 s(n1_cgal, n2_cgal);
            auto distance_squared = CGAL::squared_distance(s, v_cgal);
            d = std::sqrt(distance_squared);
        }
        else{
            auto n1 = ReferenceVertices[neighbours[0].idx-1];
            auto n2 = ReferenceVertices[neighbours[0].idx];
            auto n3 = ReferenceVertices[neighbours[0].idx + 1];
            Point_3 n1_cgal(n1[0], n1[1], n1[2]), n2_cgal(n2[0], n2[1], n2[2]), n3_cgal(n3[0], n3[1], n3[2]);
            Segment_3 s12(n1_cgal, n2_cgal), s23(n2_cgal, n3_cgal);
            auto distance_squared_12 = CGAL::squared_distance(s12, v_cgal);
            auto distance_squared_23 = CGAL::squared_distance(s23, v_cgal);
            auto distance_squared = std::min(distance_squared_12, distance_squared_23);
            d = std::sqrt(distance_squared);
        }
        accu_distance += d;
    }

    float length = 0;
    for (int i = 1; i < n; ++i) {
        auto p1 = ReferenceVertices[i - 1];
        auto p2 = ReferenceVertices[i];
        length += Vrui::Geometry::dist(p1, p2);
    }
    std::cout << "ridge::breaklineValidateProcess"
        << "the length of ReferenceBreakline is " << length
        << "\nthe accumulate distance of True positive vertices and ReferenceBreakline is " << accu_distance << std::endl;    
    //float error = accuDis / length;
}
