#include "MaPt_in_oneTrace.h"
#include <functional>
#include <queue>

using namespace masb;


inline bool MaPt_in_oneTrace::validateCandidate(float deviationAng_thres, Vector &vec1, Vector &vec2) {
    return vec1 * vec2 > deviationAng_thres;
}

size_tList MaPt_in_oneTrace::find_trace(pt_Trace_pram&power, size_t &pt_idx,
    mat_data &maData, ma_Geometry &maGeometry,intList &visitFlag){
    size_tList trace;
    trace.push_back(pt_idx);
    bool next = true;
    while (next) {
        kdtree2::KDTreeResultVector neighbours;
        Point pt = maData.ma_coords[pt_idx];
        maData.kdtree_ma_coords->r_nearest(pt, power.SearchRadius, neighbours);
        for (auto & candidate : neighbours) {
            if (visitFlag[candidate.idx] == 1) {
                next = false;
                continue;
            }
            Vector vec_pt_bis = maGeometry.ma_bisector[pt_idx];
            Point candidate_pt = maData.ma_coords[candidate.idx];
            auto vec_ptCand = candidate_pt - pt;
            vec_ptCand.normalize();
            if (maData.ma_radius[candidate.idx] < maData.ma_radius[pt_idx] &&
                validateCandidate(power.deviationAng_thres, vec_pt_bis, vec_ptCand)) {
                trace.push_back(candidate.idx);
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
void MaPt_in_oneTrace::processing(pt_Trace_pram& power, mat_data &madata, ma_Geometry &maGeometry, Sheet_idx_List &sheets) {
    //params.update();
    std::cout << "there are " << sheets.size() << " sheets\n";
    size_tList sheet_idx;
    size_t seg_id = 0;
    this->candidate_r.resize(sheets.size());
    this->candidate_cos.resize(sheets.size());
    this->candidate_dir.resize(sheets.size());
    this->segmentation.resize(sheets.size());
    this->all_traces.resize(sheets.size());
#pragma omp parallel for private(sheet_idx)
    for (size_t i = 0; i < sheets.size(); ++i) {
        sheet_idx = sheets[i];
        seg_id = i + 1;
        //std::cout << "start sheet " << seg_id << "\n";
        intList visitFlag; visitFlag.resize(sheet_idx.size(),0);
        PointList candidatept_r, candidatept_cos;
        VectorList candidatept_dir;
        traces_in_one_sheet traces_in_one_sheet;

        masb::idx_filter filter;
        mat_data maSheetData;
        ma_Geometry masheetGeometry;
        filter.processing(madata, maGeometry, sheet_idx, maSheetData, masheetGeometry);
        //intList ma_pidx; ma_pidx = sheet_idx;
        if (maSheetData.kdtree_ma_coords == NULL)
            maSheetData.kdtree_ma_coords = new kdtree2::KDTree(maSheetData.ma_coords, true);
        maSheetData.kdtree_ma_coords->sort_results = true;

        typedef std::pair<size_t, float> index_rdaius_pair;
        auto cmp = [](index_rdaius_pair left, index_rdaius_pair right) {return left.second < right.second; };
        std::priority_queue<index_rdaius_pair, std::vector<index_rdaius_pair>, decltype(cmp)> pq(cmp);
        for (size_t i = 0; i < maSheetData.ma_ptsize; ++i) {
            pq.push(index_rdaius_pair(i,maSheetData.ma_radius[i]));
        }

        while (pq.size() > 0) {
            //std::cout << "pt_idx" << pq.top().first 
            //    << "\tpt_radius " << pq.top().second << "\n";
            auto pt_idx = pq.top().first; pq.pop();
            
            if (visitFlag[pt_idx] == 0) {
                visitFlag[pt_idx] = 1;
                auto a_trace = find_trace(power, pt_idx, maSheetData, masheetGeometry, visitFlag);
                if (a_trace.size() > 1) {
                    //auto idx = a_trace.back();
                    auto pt = maSheetData.ma_coords[a_trace.back()];
                    auto r = maSheetData.ma_radius[a_trace.back()];
                    auto bisec = masheetGeometry.ma_bisector[a_trace.back()];
                    auto SepAng = masheetGeometry.ma_SeperationAng[a_trace.back()];
                    //auto vec = bisec * r;
                    //auto can = pt + bisec * r;
                    candidatept_r.push_back(pt + bisec * r);
                    //auto can = pt + bisec * (r / cos(SepAng / 2));
                    candidatept_cos.push_back(pt + bisec * (r / cos(SepAng / 2)));
                    candidatept_dir.push_back(masheetGeometry.ma_normal[a_trace.back()]);
                    
                    PointList a_trace_pt;
                    a_trace_pt.reserve(a_trace.size());
                    for (auto idx : a_trace) {
                        a_trace_pt.push_back(maSheetData.ma_coords[idx]);
                    }
                    traces_in_one_sheet.push_back(a_trace_pt);
                }
            }
        }
        //std::cout << "finish sheet " << seg_id << "\n";
        this->candidate_r[i] = candidatept_r;
        this->candidate_cos[i] = candidatept_cos;
        this->candidate_dir[i] = candidatept_dir;
        this->segmentation[i] = seg_id;
        this->all_traces[i] = traces_in_one_sheet;
        this->candidate_size += candidatept_r.size();
        
    }
}