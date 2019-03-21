#ifndef MAPT_INONETRACE_
#define MAPT_INONETRACE_
#include <madata.h>
#include <iostream>

#include "Ma_utility.h"

namespace masb{
    class pt_Trace_pram {
    public:
        std::string method = "fixe radius";
        int k_neighbours;
        float SearchRadius = 45;
        float deviationAng_thres = 10.0;
        pt_Trace_pram() {
            this->deviationAng_thres = cos((this->deviationAng_thres / 180.0)*PI);
        }
    };
    class MaPt_in_oneTrace {
    public:
        pt_Trace_pram power;
        std::vector<PointList> candidate_r, candidate_cos;
        std::vector <VectorList> candidate_dir;
        typedef std::vector<PointList> traces_in_one_sheet;
        typedef std::vector<traces_in_one_sheet> trace_shape;
        trace_shape all_traces;
        std::vector<size_t> segmentation;
        size_t candidate_size = 0;

        void processing(pt_Trace_pram& power,mat_data &madata, ma_Geometry &maGeometry, Sheet_idx_List &sheets);
    private:
        size_tList find_trace(pt_Trace_pram&power, size_t &pt_idx,
            mat_data &maData, ma_Geometry &maGeometry, intList &visitFlag);
        inline bool validateCandidate(float deviationAng_thres, Vector &vec1, Vector &vec2);
    };
}

#endif