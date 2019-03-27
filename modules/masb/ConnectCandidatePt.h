#ifndef CONNECTCANDIDATEPT_
#define CONNECTCANDIDATEPT_

#include <iostream>
#include "Ma_utility.h"

namespace ridge{
    masb::intList connectCandidatePtProcess(masb::PointList &pointcloud,masb::PointList &candidate, masb::intList &seg_id,
        masb::VectorList &direction, masb::VectorList &bisec_p, masb::VectorList &bisec_q);
}
#endif