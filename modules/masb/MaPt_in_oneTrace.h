#ifndef MAPT_INONETRACE_
#define MAPT_INONETRACE_
#include <madata.h>
#include <iostream>

namespace masb{
    struct pt_Trace_pram {
    public:
        float SearchRadius;
        float deviationAng_thres;//cos value
        pt_Trace_pram(float r = 45, float Ang_deg = 10) :SearchRadius(r), deviationAng_thres(cos((Ang_deg / 180.0)*PI)) {
            std::cout << "PtinDirectionLine_param constructor:: SearchRadius = " << SearchRadius
                << " cos(deviationAng_thres) = " << deviationAng_thres << std::endl;
        }
    };
    class MaPt_in_oneTrace {
    public:
        pt_Trace_pram power;
        ridge::PolyineList traces;
        masb::intList trace_seg_id;

        void MaPt_in_oneTrace::processing(pt_Trace_pram& power, MAT &mat, masb::intList &seg_id);
    private:
        masb::PointList find_trace(pt_Trace_pram&power, size_t &pt_idx,
            masb::PointList &cur_pt, masb::floatList &cur_radius, masb::VectorList &cur_direc, 
            kdtree2::KDTree* ptkdtree,intList &visitFlag);
        inline bool validateCandidate(float deviationAng_thres, Vector &vec1, Vector &vec2);
    };

}
#endif