#include "ConnectCandidatePt.h"


masb::intList ridge::connectCandidatePtProcess(masb::PointList &pointCloud,masb::PointList &candidate, 
    masb::intList &seg_id,masb::VectorList &direction, 
    masb::VectorList &bisec_p, masb::VectorList &bisec_q) {

    masb::intList filter;

    kdtree2::KDTree* pc_kdtree;
    pc_kdtree = new kdtree2::KDTree(pointCloud, true);
    pc_kdtree->sort_results = true;
    std::vector<float> dis_list; dis_list.reserve(candidate.size());
    for (int i = candidate.size() - 1; i >= 0; i--) {
        kdtree2::KDTreeResultVector neighbours;
        pc_kdtree->n_nearest(candidate[i], 1, neighbours);
        dis_list.push_back(neighbours[0].dis);
        if (neighbours[0].dis > 30.0) {
            filter.push_back(3);
            /*
            candidate.erase(candidate.begin() + i);
            seg_id.erase(seg_id.begin() + 1);
            direction.erase(direction.begin() + i);
            bisec_p.erase(bisec_p.begin() + i);
            bisec_q.erase(bisec_q.begin() + i);
            */
        }
        else
            filter.push_back(0);
    }
    dis_list;
    return filter;


}