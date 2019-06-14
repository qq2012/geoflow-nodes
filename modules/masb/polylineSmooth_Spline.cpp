#include "polylineSmooth.h"
#include <iostream>
#include <unsupported\Eigen\Splines>
//#include "C:\dev\Eigen3\include\eigen3\unsupported\Eigen\Splines"
//#include "C:\dev\Eigen3\include\eigen3\Eigen\Core"

#include "BSpline.h"
#include "Bezier.h"
using namespace Eigen;

ridge::PolyineList ridge::polylineSmooth8EigenSpline(ridge::PolyineList &polylines) {

    PolyineList resPolylines;
    polylines = resPolylines;
    /*
    this example return the same valus, does not interpolation
    int degree = 3;
    typedef Spline<float, 3> Spline3d;//Eigen::Spline< _Scalar, _Dim, _Degree >
    //control point size = n+1 p[0,1,2,3,...,n]
    //knot vector size = m+1   u[0,1,2,...,m]
    //Degree p = 3
    // m = n + p + 1, m = n + 4 = control point size + 3
    typedef Spline3d::PointType PointType;
    typedef Spline3d::KnotVectorType KnotVectorType;
    typedef Spline3d::ControlPointVectorType ControlPointVectorType;

    for (auto &apolyline : polylines) {
        //use polyline points as control points, and uniform spline

        int pt_size = apolyline.size();
        Eigen::MatrixXf matrix_pt(3, pt_size);
        int col = 0;
        for (auto& p : apolyline) {
            matrix_pt.col(col++) << p[0], p[1], p[2];
            //col++;
        }
        //std::cout << matrix_pt.col(0);
        //std::cout << "\nThe matrix matrix_pt is of size "
        //    << matrix_pt.rows() << "x" << matrix_pt.cols() << std::endl;
        //std::cout << "It has " << matrix_pt.size() << " coefficients" << std::endl;

        ControlPointVectorType ctrls_points;// = ControlPointVectorType::Random(3, 100);//Random(hang,lie)
        ControlPointVectorType points = ControlPointVectorType::Random(3, 20);
        //
        //ctrls_points = [ x0,x1,x2,...,xn
        //                 y0,y1,y2,...,yn
        //                 z0,z1,z2,...,zn ]
        //ctrls_points = matrix_pt;
        ctrls_points = points;
        for (int i = 0; i < ctrls_points.cols(); ++i) {
            PointType pt = ctrls_points.col(i);
            testpolyline.push_back(masb::Point(pt(0), pt(1), pt(2)));
        }

        const Spline3d spline = SplineFitting<Spline3d>::Interpolate(ctrls_points, degree);
        //Interpolate(const PointArrayType & 	pts,
        //                  DenseIndex 	degree,
        //    const KnotVectorType & 	knot_parameters)

        KnotVectorType chord_lengths; // knot parameters
        Eigen::ChordLengths(ctrls_points, chord_lengths);

        masb::PointList cur_ans_line;
        std::cout << "ctrls_points.cols()" << ctrls_points.cols() << std::endl;
        for (Eigen::DenseIndex i = 0; i < ctrls_points.cols(); ++i) {
            PointType pt = spline(chord_lengths(i));
            PointType ref = ctrls_points.col(i);
            if ((pt - ref).matrix().norm() < 1e-14)
                std::cout << "(smoothed pt - raw pt).norm()< 1e-14 at idx -- " << i << std::endl;
            //VERIFY((pt - ref).matrix().norm() < 1e-14);
            //VERIFY() serves the same purpose as ASSERT()
            std::cout << pt.col(0);
            std::cout << "\npt(0)" << pt(0) << ",pt(1)" << pt(1) << ",pt(2)" << pt(2);
            std::cout << "\npt[0]" << pt[0] << ",pt[1]" << pt(1) << ",pt[2]" << pt[2];
            std::cout << "\nThe point matrix pt is of size "
                << pt.rows() << "x" << pt.cols() << std::endl;
            std::cout << "It has " << pt.size() << " coefficients" << std::endl;
            //auto x = pt(0); matrix(i)--i raw (hang)
            //auto t = pt[0];
            cur_ans_line.push_back(masb::Point(pt(0), pt(1), pt(2)));
        }
        resPolylines.push_back(cur_ans_line);

    */

    /*
    this example return infinite interpolation value, but suessfully return smooth line for symple line
    for (auto &aline : polylines) {
        auto z = aline[0][2];
        masb::PointList a_resPolylines;
        int size = aline.size();
        Eigen::VectorXd xvals(size);
        Eigen::VectorXd yvals(size);
        for (int i = 0; i < size; ++i) {
            auto pt = aline[i];
            xvals(i) = pt[0];
            yvals(i) = pt[1];
            //std::cout <<"pt[0] "<< pt[0]<< ", xvals(0) " << xvals(i)
            //    << "pt[1]"<<pt[1]<<"yvals(0) " << yvals(i) << std::endl;
        }
        for (int i = 0; i < size; ++i) {
            
        }
        ridge::SplineFunction s(xvals, yvals);
        int step_num = size;
        float step = (s.x_max - s.x_min) / step_num;
        for (int i = 0; i < step_num-1; i++) {
        //for (float x = s.x_min; x < s.x_max; x+=10) {
            float x = s.x_min + step * i;
            auto y = (float)s(x);
            a_resPolylines.push_back(masb::Point(x, y, z));
            std::cout << "x--" << x << ",y--" << s(x) << ",z " << z << std::endl;
        }
        resPolylines.push_back(a_resPolylines);
    }
    */
    masb::PointList a_resPolylines, testpolyline;

    Eigen::VectorXd xvals(3);
    Eigen::VectorXd yvals(xvals.rows());

    xvals << 0, 15, 30;
    yvals << 0, 12, 17;
    //auto a = xvals(1);
    std::cout << "xvals(0) " << xvals(0)
        << ",xvals(1) " << xvals(1)
        << ",xvals(2) " << xvals(2) << "\n" << std::endl;

    testpolyline.push_back(masb::Point((float)xvals(0), (float)yvals(0), 0));
    testpolyline.push_back(masb::Point((float)xvals(1), (float)yvals(1), 0));
    testpolyline.push_back(masb::Point((float)xvals(2), (float)yvals(2), 0));
    polylines.push_back(testpolyline);

    ridge::SplineFunctionEigen s(xvals, yvals);
    int step_num = 10;
    float step = (s.x_max - s.x_min) / step_num;
    for (int i = 0; i < step_num; i++) {
        float x = s.x_min + step * i;
        auto y = (float)s(x);
        a_resPolylines.push_back(masb::Point(x, y, 0));
        std::cout << s(x) << std::endl;
    }
    resPolylines.push_back(a_resPolylines);
    
    return resPolylines;
    
}

ridge::PolyineList ridge::polylineBSplineSmooth(ridge::PolyineList &polylines) {
    std::cout << "\nridge::polylineBSplineSmooth()" << std::endl;
    PolyineList smoothedPolylines;
    for (auto &a_polyline : polylines) {

        std::cout << "before smooth a_polyline.size = " << a_polyline.size() << std::endl;

        int size = a_polyline.size()*1.5;
        masb::PointList aSmothLine;
        Curve* curve = new Bezier();
        curve->set_steps(size); // generate size interpolate points fot the last way points
        for (auto& pt : a_polyline) {
            curve->add_way_point(Vector(pt[0], pt[1], pt[2]));
        }
        //std::cout << "nodes: " << curve->node_count() << std::endl;
        //std::cout << "total length: " << curve->total_length() << std::endl;
        for (int i = 0; i < curve->node_count(); ++i) {
            //std::cout << "node #" << i << ": " << curve->node(i).toString()
            //    << " (length so far: " << curve->length_from_starting_point(i) << ")" << std::endl;
            auto pt = curve->node(i);
            aSmothLine.push_back(masb::Point((float)pt.x, (float)pt.y, (float)pt.z));
        }
        smoothedPolylines.push_back(aSmothLine);
        delete curve;

        std::cout << "after smooth aSmothLine.size = " << aSmothLine.size() << std::endl;
    }
    return smoothedPolylines;
    /*
    int size = polylines[0].size()*1.5;
    PolyineList resPolylines;
    polylines = resPolylines;

    Curve* curve = new Bezier();
	curve->set_steps(size); // generate 100 interpolate points between the last 4 way points
    curve->add_way_point(Vector(55.5, 37, 2));
    curve->add_way_point(Vector(43, 1.8, 22));
    curve->add_way_point(Vector(1.3, 1.4, 10));
    curve->add_way_point(Vector(24.5, 73, 42));
    curve->add_way_point(Vector(38, 28, 50));
	curve->add_way_point(Vector(43, 61.8, 30));

    masb::PointList a_resPolylines, testpolyline;
    testpolyline.push_back(masb::Point(1.3, 1.4, 10));
    testpolyline.push_back(masb::Point(24.5, 73, 42));
    testpolyline.push_back(masb::Point(38, 28, 50));
    testpolyline.push_back(masb::Point(43, 61.8, 30));
    polylines.push_back(testpolyline);

	std::cout << "nodes: " << curve->node_count() << std::endl;
	std::cout << "total length: " << curve->total_length() << std::endl;
	for (int i = 0; i < curve->node_count(); ++i) {
		std::cout << "node #" << i << ": " << curve->node(i).toString() 
            << " (length so far: " << curve->length_from_starting_point(i) << ")" << std::endl;
        auto pt = curve->node(i);
        a_resPolylines.push_back(masb::Point((float)pt.x, (float)pt.y, (float)pt.z));
	}
    resPolylines.push_back(a_resPolylines);
    */
	
    //return resPolylines;
}