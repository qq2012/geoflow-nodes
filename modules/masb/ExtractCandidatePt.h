#ifndef EXTRACTCANDIDATEPT_
#define EXTRACTCANDIDATEPT_

#include <iostream>
#include "madata.h"
//using namespace masb;
namespace masb {
    struct ExtractCandidatePt_pram {
        std::string method;

    };

    class ExtractCandidatePt {
    public:
        int i;
        void processing();
    };

    class CandidateDir {
    public:
        VectorList CandidateDirection;
        void estimateCandidateDir(ma_data &madata, ma_Geometry &maGeometry);
    };
}
#endif