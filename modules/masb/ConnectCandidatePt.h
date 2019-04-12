#ifndef CONNECTCANDIDATEPT_
#define CONNECTCANDIDATEPT_

#include <iostream>
#include "Ma_utility.h"

namespace ridge{
    typedef std::pair<masb::Point, masb::Point> PointPair;
    typedef std::vector<PointPair> segment;
    typedef std::vector<masb::PointList> line;
    void connectCandidatePt8MST(masb::PointList &PointCloud,masb::PointList &candidate, 
        masb::intList &seg_id, masb::intList &filter, segment &segmentList, masb::intList &idList,
        line &symple_segmentList, masb::intList &symple_idList);
    void connectCandidatePt8MST_nosegid(masb::PointList &PointCloud, masb::PointList &candidate, masb::intList &filter, segment &segmentList);
    void connectCandidatePtSmooth(line &symple_segmentList, line &smoothLine);
    void FindTopology(line &smoothLines);
    void connectCandidatePtSmooth8Polynomials(line &Lines);
    void connectCandidatePt8Spline(masb::PointList &PointCloud, masb::PointList &candidate,
        masb::intList &seg_id, masb::intList &filter, ridge::line &line_segmentList, masb::intList &idList);
}
#endif