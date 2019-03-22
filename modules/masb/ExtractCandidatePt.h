#ifndef EXTRACTCANDIDATEPT_
#define EXTRACTCANDIDATEPT_

#include <iostream>
#include "madata.h"
#include "Ma_utility.h"

//using namespace masb;
namespace masb {
    /*
    struct ExtractCandidatePt_pram {
        std::string method;

    };
    
    class ExtractCandidatePt {
    public:
        int i;
        void processing();
    };
    */
    class ExtractCandidatePt_pram {
    public:
        std::string method = "fixe radius";
        int k_neighbours;
        float SearchRadius = 45;
        float deviationAng_thres = 10.0;
        ExtractCandidatePt_pram() {
            this->deviationAng_thres = cos((this->deviationAng_thres / 180.0)*PI);
        }
    };
    class ExtractCandidatePt {
    public:
        ExtractCandidatePt_pram power;
        std::vector<PointList> candidate_r, candidate_cos;
        std::vector <VectorList> candidate_dir;
        size_t candidate_size = 0;

        void processing(ExtractCandidatePt_pram & power, mat_data &madata, ma_Geometry &maGeometry, Sheet_idx_List &sheets);
    private:
        size_tList find_trace(ExtractCandidatePt_pram &power, size_t &pt_idx,
            mat_data &maData, ma_Geometry &maGeometry, intList &visitFlag);
        inline bool validateCandidate(float deviationAng_thres, Vector &vec1, Vector &vec2);
    };
}
#endif