#ifndef CONNECTCANDIDATEPT_
#define CONNECTCANDIDATEPT_

#include <iostream>
#include "Ma_utility.h"

namespace ridge{
    typedef std::pair<masb::Point, masb::Point> PointPair;
    typedef std::vector<PointPair> segment;
    typedef std::vector<masb::PointList> line;
    
    typedef std::pair<int, int> int_pair;
    typedef std::vector<int_pair> int_pair_vec;

    void connectCandidatePt8MST(masb::PointList &PointCloud,masb::PointList &candidate, 
        masb::intList &seg_id, masb::intList &filter, segment &segmentList, masb::intList &idList,
        line &symple_segmentList, masb::intList &symple_idList);
    void connectCandidatePt8MST_nosegid(masb::PointList &PointCloud, masb::PointList &candidate, masb::intList &filter, segment &segmentList);
    void connectCandidatePtSmooth(line &symple_segmentList, line &smoothLine);
    void FindTopology(line &smoothLines, const masb::intList &symple_idList, const int_pair_vec &adjacency);
    void connectCandidatePtSmooth8Polynomials(line &Lines);
    void connectCandidatePt8Spline(masb::PointList &PointCloud, masb::PointList &candidate,
        masb::intList &seg_id, masb::intList &filter, ridge::line &line_segmentList, masb::intList &idList);
}
#endif