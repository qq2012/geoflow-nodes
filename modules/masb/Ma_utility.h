#ifndef MA_UTILITY_
#define MA_UTILITY_

#include "madata.h"
namespace masb {

    class idx_filter {
    public:
        void processing(ma_data &madata, ma_Geometry &maGeometry,intList &remaining_idx, 
            mat_data &remainingData, intList &remainingma_in_out, ma_Geometry &remainingGeometry);
        void processing(ma_data &madata, intList &remaining_idx,
            mat_data &remainingData, intList &remainingma_in_out);
    };

    class MaPt_in_oneTrace {
        typedef std::vector<intList> pt_Traces;
    public:
        float deviationAng_thres = 10.0;
        pt_Traces traces;
        void processing(PointList &ma_coords, ma_Geometry &maGeometry, pt_Traces &traces);
        MaPt_in_oneTrace() {
            this->deviationAng_thres = cos((this->deviationAng_thres / 180.0)*PI);
        }
    private:
        float caldeviationAng(Point &pt1, Vector &pt1_norm, Point &pt2);
    };

    
}

#endif