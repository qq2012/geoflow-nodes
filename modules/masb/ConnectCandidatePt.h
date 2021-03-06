#ifndef CONNECTCANDIDATEPT_
#define CONNECTCANDIDATEPT_

#include <iostream>
#include "Ma_utility.h"

namespace ridge{
    /*
    typedef std::pair<masb::Point, masb::Point> PointPair;
    typedef std::vector<PointPair> LineSegmentList;
    typedef std::vector<masb::PointList> PolyineList;
    
    typedef std::pair<int, int> int_pair;
    typedef std::vector<int_pair> int_pair_vec;
    */
    void connectCandidatePt8MST(masb::PointList &pts, masb::VectorList &pt_directon, masb::intList &pt_id,
        LineSegmentList &mstLineSegment, masb::intList &mstLineSegment_id, 
        ridge::PolyineList &polylines_maxDistance, ridge::PolyineList &polylines_maxAccDist, ridge::PolyineList &polylines_maxPts, 
        masb::intList &polyline_id);
    //void connectCandidatePt8MST_nosegid(masb::PointList &PointCloud, masb::PointList &candidate, masb::intList &filter, segment &segmentList);
    void FindTopology(PolyineList &polylines, const masb::intList &polylinesID, const int_pair_vec &adjacency);

    void ConnectCandidatePt8PolynomialFitting(masb::PointList &pts, masb::intList &pt_id, masb::PointList &pointCloud, 
        float error_thresh, ridge::PolyineList &polylines, const int &order);

    }
#endif