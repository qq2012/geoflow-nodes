#include "Ma_utility.h"
#include <iostream>

//#include <Eigen/Eigenvalues>

// Vrui
#include <vrui/Geometry/ComponentArray.h>
#include <vrui/Math/Math.h>
#include <vrui/Geometry/PCACalculator.h>

#include <Eigen/Geometry>

void masb::coordinateTransformation_2d::ComputeNewAxis(const masb::PointList &input_pts) {
    /*
    // create matrix
    Eigen::MatrixXd A(input_pts.size(), 2);
    for (size_t i = 0; i < input_pts.size(); ++i) {
        auto p = input_pts[i];
        A(i, 0) = p[0];
        A(i, 1) = p[1];
    }
    //A=[[x1,y1],[x2,y2],[x3,y3].......]

    Eigen::EigenSolver<Eigen::MatrixXd> es(A);
    std::cout << "The eigenvalues of A are:" << std::endl << es.eigenvalues() << std::endl;
    std::cout << "The matrix of eigenvectors, V, is:" << std::endl << es.eigenvectors() << std::endl << std::endl;

    std::complex<double> lambda = es.eigenvalues()[0];
    std::cout << "Consider the first eigenvalue, lambda = " << lambda << std::endl;
    Eigen::VectorXcd v = es.eigenvectors().col(0);
    std::cout << "If v is the corresponding eigenvector, then lambda * v = " << std::endl << lambda * v << std::endl;
    std::cout << "... and A * v = " << std::endl << A.cast<std::complex<double> >() * v << std::endl << std::endl;
    Eigen::MatrixXcd D = es.eigenvalues().asDiagonal();
    Eigen::MatrixXcd V = es.eigenvectors();
    std::cout << "Finally, V * D * V^(-1) = " << std::endl << V * D * V.inverse() << std::endl;
    */

    Vrui::Geometry::PCACalculator<2> PCACalc;

    for (auto&p:input_pts)
        PCACalc.accumulatePoint(Point_2d( p[0],p[1] ));

    double eigen_values[2];
    PCACalc.calcCovariance();
    PCACalc.calcEigenvalues(eigen_values);
    this->x_axis = PCACalc.calcEigenvector(eigen_values[0]);
    this->y_axis = PCACalc.calcEigenvector(eigen_values[1]);

}

void masb::coordinateTransformation_2d::Transform(const masb::PointList &input_pts, masb::PointList &Transform_Pts) {

    masb::coordinateTransformation_2d::ComputeNewAxis(input_pts);

    Eigen::MatrixXd TransformMatrix(2, 2);
    TransformMatrix(0, 0) = this->x_axis[0];
    TransformMatrix(0, 1) = this->x_axis[1];
    TransformMatrix(1, 0) = this->y_axis[0];
    TransformMatrix(1, 1) = this->y_axis[1];

    //Eigen::Transform<float, 2, Eigen::Affine> t;
    //t.data = TransformMatrix;
    Transform_Pts.resize(input_pts.size());

    for (auto &p : input_pts) {
        Eigen::MatrixXd p_xy(2, 1);
        p_xy(0, 0) = p[0];
        p_xy(1, 0) = p[1];
        auto p_t = TransformMatrix * p_xy;
        Transform_Pts.push_back(masb::Point(p_t(0, 0), p_t(1, 0), p[2]));
    }

}

void masb::coordinateTransformation_2d::InverseTransform(const masb::PointList &pts, masb::PointList &Transform_Pts) {

    Eigen::MatrixXd TransformMatrix(2, 2);
    TransformMatrix(0, 0) = this->x_axis[0];
    TransformMatrix(0, 1) = this->x_axis[1];
    TransformMatrix(1, 0) = this->y_axis[0];
    TransformMatrix(1, 1) = this->y_axis[1];

    Transform_Pts.resize(pts.size());

    for (auto &p : pts) {
        Eigen::MatrixXd p_xy(2, 1);
        p_xy(0, 0) = p[0];
        p_xy(1, 0) = p[1];
        auto p_t = TransformMatrix.reverse() * p_xy;
        Transform_Pts.push_back(masb::Point(p_t(0, 0), p_t(1, 0), p[2]));
    }
}
