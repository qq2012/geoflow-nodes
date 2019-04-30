#include "breaklineValidate.h"

#include <CGAL/Simple_cartesian.h>
typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_3 Point_3;
typedef Kernel::Segment_3 Segment_3;

float ridge::breaklineValidateProcess(ridge::line &TestLine, ridge::line &TruthLine) {
    Point_3 p(1.1, 2.2, 3.3);
    Point_3 n1(0, 0, 0), n2(1, 1, 1);
    Segment_3 s(n1, n2);
    auto distance = CGAL::squared_distance(s, p);
    float accuDis ,length;
    float error = accuDis / length;
    return error;
}
