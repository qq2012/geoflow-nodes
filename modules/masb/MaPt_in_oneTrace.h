#ifndef MAPT_INONETRACE_
#define MAPT_INONETRACE_
#include <madata.h>
#include <iostream>
namespace masb{
    struct pt_in_oneTrace_pram {
        std::string method = "fixe radius";
        int k_neighbours;
        float searchRadius;
        float deviationAng_thres = 10.0;
        pt_in_oneTrace_pram() {
            this->deviationAng_thres = cos((this->deviationAng_thres / 180.0)*PI);
        }
    };
    class MaPt_in_oneTrace {
        typedef std::vector<intList> pt_Traces;
    public:
        pt_in_oneTrace_pram params;
        pt_Traces traces;
        void processing(PointList &ma_coords, ma_Geometry &maGeometry, pt_Traces &traces);
    private:
        float caldeviationAng(Point &pt1, Vector &pt1_norm, Point &pt2);
    };
}

#endif