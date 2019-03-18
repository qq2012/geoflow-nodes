#ifndef MA_UTILITY_
#define MA_UTILITY_

#include "madata.h"
#include <variant>
namespace masb {

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
}

#endif