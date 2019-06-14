#ifndef PTINDIRECTIONLine_
#define PTINDIRECTIONLine_

#include "madata.h"
#include "Ma_utility.h"

namespace ridge {
    struct PtinDirectionLine_param{
    public:
        float SearchRadius;
        float deviationAng_thres;//cos value
        PtinDirectionLine_param(float r=45, float Ang_deg=10):SearchRadius(r),deviationAng_thres(cos((Ang_deg / 180.0)*PI)){
            std::cout << "PtinDirectionLine_param constructor:: SearchRadius = " << SearchRadius
                << " cos(deviationAng_thres) = " << deviationAng_thres << std::endl;
        }
    };

    class PtinDirectionLine {
    public:
        ridge::PolyineList lines;
        masb::intList lines_id;
        void processing(PtinDirectionLine_param& power, masb::PointList &pts, masb::intList& pt_id, masb::VectorList &vectors);
    private:
        inline bool validateCandidate(float deviationAng_thres, masb::Vector &vec1, masb::Vector &vec2);
    };
}
#endif // !PTINDIRECTIONLine_
