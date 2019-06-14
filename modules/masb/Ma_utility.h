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

    class idx_filter {
    public:
        void processing(ma_data &madata, ma_Geometry &maGeometry,intList &remaining_idx, 
            mat_data &remainingData, ma_Geometry &remainingGeometry);
        void processing(mat_data &madata, ma_Geometry &maGeometry, size_tList &remaining_idx,
            mat_data &remainingData, ma_Geometry &remainingGeometry);
        void processing(ma_data &madata, intList &remaining_idx,mat_data &remainingData);
    };

    class splitInOut {
    private:
        size_t size ;
        size_t begin ;
    public:
        void processing(mat_data &madata, ma_Geometry &maGeometry, bool InFlag,
            mat_data &madata_new, ma_Geometry &maGeometry_new);
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