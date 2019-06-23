#ifndef MA_UTILITY_
#define MA_UTILITY_

#include "madata.h"
#include <variant>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>

#include <sstream>
#include <algorithm>
#include <iterator>
#include <vector>
#include <iterator>

namespace masb {

    void PolynomialFitting(const masb::PointList &pts, const int &order, masb::PointList &polyline,float &min_error);

    struct condition_pram {
    public:
        float MaxEdgeBallRadius; //filter edge ball with big r
        float MinEdgeBallRadius; //filter edge ball with small r
        float filterDistance; //filter edge ball far away from ground
        float unshrinkingDist;//distance from candidate pt (x,y,z_interpolation) to unshrinking point 
        condition_pram(float r_x, float r_n, float d_p, float d_u) :MaxEdgeBallRadius(r_x),
            MinEdgeBallRadius(r_n), filterDistance(d_p), unshrinkingDist(d_u) {
            std::cout << "condition_pram constructor" << std::endl;
        }
    };
    intList conditionFilter(masb::MAT &mat, masb::condition_pram& power, PointList &pointcloud,
        PointList &unShrinkingPt, intList &seg_id);

    class coordinateTransformation_2d {
    public:
        masb::Vector_2d x_axis;
        masb::Vector_2d y_axis;
        void Transform(const masb::PointList &input_pts, masb::PointList &Transform_Pts);
        void InverseTransform(const masb::PointList &pts, masb::PointList &Transform_Pts);
    private:
        void ComputeNewAxis(const masb::PointList &input_pts);
    };

    //typedef std::variant<int, float, Point, Vector> vectype;
    template <typename T>
    inline std::vector<T> VectorSlice(std::vector<T> & oldVec,
        size_t const size, size_t const start) {
        typedef std::vector<T>::const_iterator iterator_typ;
        iterator_typ first = oldVec.begin() + start;
        iterator_typ last = oldVec.begin() + start + size;
        std::vector<T> newVec(first, last);
        return newVec;
    }
    template<typename Out>
    void split(const std::string &s, char delim, Out result) {
        std::stringstream ss(s);
        std::string item;
        while (std::getline(ss, item, delim)) {
            *(result++) = item;
        }
    }
    inline std::vector<std::string> split2nums(const std::string &s, char delim) {
        std::vector<std::string> elems;
        split(s, delim, std::back_inserter(elems));
        return elems;
    }


}
#endif