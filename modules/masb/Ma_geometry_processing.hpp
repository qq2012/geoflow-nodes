#ifndef MA_GEOMETRY_PROCESSING_
#define MA_GEOMETRY_PROCESSING_

#include "madata.h"
#include "Ma_utility.h"

namespace masb {
    struct maGeom_result {
        Vector bisector;
        float SepAng;
        Vector norm;
    };
    void compute_ma_geometry(ma_data &madata, ma_Geometry &maGeometry);
    //void compute_ma_geometry(ma_data &madata, ma_Geometry &maGeometry, intList &remaining_idx);
    void compute_ma_geometry_new(const ma_data &madata, MAT &mat,const bool isInterior);
}
#endif // !MA_GEOMETRY_PROCESSING_
