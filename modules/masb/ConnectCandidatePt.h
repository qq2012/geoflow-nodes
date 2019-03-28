#ifndef CONNECTCANDIDATEPT_
#define CONNECTCANDIDATEPT_

#include <iostream>
#include "Ma_utility.h"

namespace ridge{
    typedef std::pair<masb::Point, masb::Point> PointPair;
    typedef std::vector<PointPair> segment;
    void connectCandidatePt8MST(masb::PointList &PointCloud,masb::PointList &candidate, 
        masb::intList &seg_id, masb::intList &filter, segment &segmentList, masb::intList &idList);
}
#endif