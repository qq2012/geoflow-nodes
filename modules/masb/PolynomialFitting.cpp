#include "Ma_utility.h"

#include <unsupported/Eigen/Polynomials>

void polyfit(const std::vector<double> &xv, const std::vector<double> &yv, std::vector<double> &coeff,
    float &min_error,int order){

    Eigen::MatrixXd A(xv.size(), order + 1);
    //auto t = yv.front(); //double
    //yv.size() how many element in yv;
    Eigen::VectorXd yv_mapped = Eigen::VectorXd::Map(&yv.front(), yv.size());
    Eigen::VectorXd result;

    assert(xv.size() == yv.size());
    assert(xv.size() >= order + 1);

    // create matrix
    for (size_t i = 0; i < xv.size(); i++)
        for (size_t j = 0; j < order + 1; j++)
            A(i, j) = pow(xv.at(i), j);

    // solve for linear least squares fit
    result = A.householderQr().solve(yv_mapped);
    
    coeff.resize(order + 1);
    for (size_t i = 0; i < order + 1; i++)
        coeff[i] = result[i];
     
    auto y_e =  A * result;
    auto ev = y_e - yv_mapped;
    auto e2 = (ev.transpose())*ev;
    /*
    std::cout << "e2.size()" << e2.size() << " e2(0, 0) = " << e2(0, 0)
        << " sqrt(e2(0, 0)) = " << sqrt(e2(0, 0))
        << " pt.size = " << xv.size()
        << " sqrt(e2(0, 0))/ pt.size = " << sqrt(e2(0, 0)) / xv.size()
        << std::endl;
        */
    min_error = sqrt(e2(0, 0));
}

void masb::PolynomialFitting(const masb::PointList &pts, masb::PointList &polyline, float &min_error) {
    
    std::vector<double> x_values, y_values, coeff;
    //float min_error;

    for (auto &pt : pts) {
        x_values.push_back(pt[0]);
        y_values.push_back(pt[1]);
    }
    
    polyfit(x_values, y_values, coeff, min_error,3);
    //printf("%f + %f*x + %f*x^2 + %f*x^3\n", coeff[0], coeff[1], coeff[2], coeff[3]);
    auto x_min_idx = std::min_element(x_values.begin(), x_values.end()) - x_values.begin();
    auto x_min = x_values[x_min_idx];
    auto x_max_idx = std::max_element(x_values.begin(), x_values.end()) - x_values.begin();
    auto x_max = x_values[x_max_idx];
    for (float x_n = x_min; x_n < x_max; ++x_n) {
        auto y_n = coeff[0] + coeff[1] * x_n + coeff[2] * pow(x_n, 2) + coeff[3] * pow(x_n, 3);
        polyline.push_back(masb::Point(x_n, y_n, 0.0));
    }
    auto y_n = coeff[0] + coeff[1] * x_max + coeff[2] * pow(x_max, 2) + coeff[3] * pow(x_max, 3);
    polyline.push_back(masb::Point(x_max, y_n, 0.0));

    //return min_error;
}