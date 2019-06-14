#include "MaPt_in_oneTrace.h"
#include <functional>
#include <queue>
#include <map>
using namespace masb;


inline bool MaPt_in_oneTrace::validateCandidate(float deviationAng_thres, Vector &vec1, Vector &vec2) {
    return vec1 * vec2 > deviationAng_thres;
}

masb::PointList MaPt_in_oneTrace::find_trace(pt_Trace_pram&power, size_t &pt_idx,
    masb::PointList &cur_pt, masb::floatList &cur_radius, masb::VectorList &cur_direc, 
    kdtree2::KDTree* ptkdtree,intList &visitFlag){

    masb::PointList trace;
    
    trace.push_back(cur_pt[pt_idx]);
    bool next = true;
    while (next) {
        kdtree2::KDTreeResultVector neighbours;
        Point pt = cur_pt[pt_idx];
        ptkdtree->r_nearest(pt, power.SearchRadius, neighbours);

        for (auto & candidate : neighbours) {
            if (visitFlag[candidate.idx] == 1) {
                next = false;
                continue;
            }
            Vector vec_pt_bis = cur_direc[pt_idx];
            Point candidate_pt = cur_pt[candidate.idx];
            auto vec_ptCand = candidate_pt - pt;
            vec_ptCand.normalize();

            if (cur_radius[candidate.idx] < cur_radius[pt_idx] &&
                validateCandidate(power.deviationAng_thres, vec_pt_bis, vec_ptCand)) {

                trace.push_back(candidate_pt);
                pt_idx = candidate.idx;
                visitFlag[candidate.idx] = 1;
                next = true;
                break;
            }
            next = false;
        }
    }
    return trace;
}
void MaPt_in_oneTrace::processing(pt_Trace_pram& power, MAT &mat, masb::intList &seg_id) {
    
    std::map<int, int> seg_frequency;
    for (auto i : seg_id) {
        ++seg_frequency[i];
    }
    /*
    for (const auto& e : seg_frequency)
        std::cout << "Element " << e.first
         << " encountered " << e.second << " times\n";
    */
    for (const auto& e : seg_frequency) {
        auto cur_sheet = e.first;
        auto cur_size = e.second;
        if (cur_sheet == 0 || cur_size < 3)
            continue;

        std::cout << "cur_sheet ID = " << cur_sheet
            << "cur_size " << e.second << std::endl;

        masb::PointList cur_pt;
        masb::floatList cur_radius;
        masb::VectorList cur_direc;
        cur_pt.reserve(cur_size);
        cur_radius.reserve(cur_size);
        cur_direc.reserve(cur_size);
        for (int i = 0; i < seg_id.size(); i++) {
            if (seg_id[i] == cur_sheet) {
                cur_pt.push_back(mat.atom[i]);
                cur_radius.push_back(mat.radius[i]);
                cur_direc.push_back(mat.bisector[i]);
            }
        }
        if (cur_pt.size() != cur_size)
            std::cout << "Error" << std::endl;

        intList cur_visitFlag; 
        cur_visitFlag.resize(cur_size,0);
        
        auto cur_pt_kdtree = new kdtree2::KDTree(cur_pt, true);
        cur_pt_kdtree->sort_results = true;

        typedef std::pair<size_t, float> index_rdaius_pair;
        auto cmp = [](index_rdaius_pair left, index_rdaius_pair right) {return left.second < right.second; };
        std::priority_queue<index_rdaius_pair, std::vector<index_rdaius_pair>, decltype(cmp)> pq(cmp);

        for (size_t i = 0; i < cur_size; ++i) {
            pq.push(index_rdaius_pair(i, cur_radius[i]));
        }

        while (pq.size() > 0) {
            //std::cout << "pt_idx" << pq.top().first 
            //    << "\tpt_radius " << pq.top().second << "\n";
            auto pt_idx = pq.top().first; 
            pq.pop();
            
            if (cur_visitFlag[pt_idx] == 0) {
                cur_visitFlag[pt_idx] = 1;
                auto a_trace = find_trace(power, pt_idx, cur_pt, cur_radius, cur_direc, 
                    cur_pt_kdtree, cur_visitFlag);
                if (a_trace.size() > 1) {
                    this->traces.push_back(a_trace);
                    this->trace_seg_id.push_back(cur_sheet);
                }
            }
        }        
    }
}