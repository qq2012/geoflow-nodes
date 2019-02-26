#include "Ma_Segmentation_processing.h"
#include <iostream>
#include <cmath>
#include <stack>
using namespace masb;
using namespace std;

//static const double PI = 3.14159265358979323846264338327950288;

void MaSeg_power::update() {
    this->bisec_thres = cos(this->bisecavg_thres / 180.0)*PI;
    this->bisecavg_thres = cos(this->bisecavg_thres / 180.0)*PI;
    this->bisecdiff_thres = cos(this->bisecdiff_thres / 180.0)*PI;
    this->theta_thres = this->theta_thres / 180.0 *PI;
    this->spokecross_thres = cos(this->spokecross_thres / 180.0)*PI;
}
inline bool MaSegProcess::if_all_segmented() {
    this->point_segment_idx;
    if (std::find(point_segment_idx.begin(), point_segment_idx.end(), -1) != point_segment_idx.end())
        return false;/* point_segment_idx -- contains -1 */
    else
        return true;/* point_segment_idx -- does not contain -1 */
}
inline size_t MaSegProcess::findseed8r(floatList *ma_radius) {
    for (size_t i = 1; i < this->size; ++i) {
        if (point_segment_idx[i] != -1 && (*ma_radius)[i] < power.seed_radius_thres)
            return i;
    }
}
inline size_t MaSegProcess::findseed() {
    for (size_t i = 1; i < this->size; ++i) {
        if (point_segment_idx[i] != -1)
            return i;
    }
}
inline bool MaSegProcess::valid_candidate_bisec(size_t idx1, size_t idx2, ma_Geometry &maGeometry) {
    Vector bisec1 = (*maGeometry.ma_bisector)[idx1];
    Vector bisec2 = (*maGeometry.ma_bisector)[idx2];
    return bisec1 * bisec2 < power.bisec_thres;
}

bool MaSegProcess::validateCandidate(size_t idx1, size_t idx2, 
    ma_data &madata, ma_Geometry &maGeometry) {
    if (power.method == "bisec")
        return valid_candidate_bisec(idx1,idx2, maGeometry);
    else
        std::cout << "not having this segmentation method yet, plz check.\n";
}

void MaSegProcess::grow_sheet(size_t initial_seed_idx, ma_data &madata, ma_Geometry &maGeometry) {
    stack<size_t>seeds;
    seeds.push(initial_seed_idx);
    vector<size_t> idx_in_sheet;
    while (!seeds.empty()) {
        size_t seed_idx = seeds.top();// .first;???????
        seeds.pop();
        Point seed_pt = (*madata.ma_coords)[seed_idx];
        kdtree2::KDTreeResultVector neighbours;
        madata.kdtree_ma_coords->n_nearest(seed_pt, power.k_neib, neighbours);
        for (auto &candidate : neighbours) {
            if (point_segment_idx[candidate.idx] != -1 && 
                validateCandidate(seed_idx, candidate.idx, madata, maGeometry)) {
                idx_in_sheet.push_back(candidate.idx);
                seeds.push(candidate.idx);
                point_segment_idx[candidate.idx] = sheet_counter;
            }
        }
    }
    if (idx_in_sheet.size() < power.mincount) {
        point_segment_idx[initial_seed_idx] = 0;
        for (size_t &idx : idx_in_sheet) {
            point_segment_idx[idx] = -1;
        }
    }
}

void MaSegProcess::processing(ma_data &madata, intList &remainingma_in_out,ma_Geometry &maGeometry) {
    this->power.update();
    this->size = madata.m * 2;
    point_segment_idx.resize(size, -1);
    
    if (madata.kdtree_ma_coords == NULL)
        madata.kdtree_ma_coords = new kdtree2::KDTree((*madata.ma_coords), true);
    madata.kdtree_ma_coords->sort_results = true;
    size_t initial_seed_idx = 0;
    while (initial_seed_idx) {
        grow_sheet(initial_seed_idx,madata, maGeometry);
        sheet_counter++;
        if (!if_all_segmented())
            initial_seed_idx = findseed();
        else
            initial_seed_idx = NULL;
    }
    std::cout << sheet_counter << '\n';
}
