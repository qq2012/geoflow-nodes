#include "Ma_Segmentation_processing.h"
#include <iostream>
#include <cmath>
#include <stack>
#include <numeric>  
#include <random>
#include <algorithm> 

using namespace masb;
using namespace std;

inline bool MaSegProcess::if_all_segmented() {
    if (std::find(point_segment_idx.begin(), point_segment_idx.end(), -1) != point_segment_idx.end())
        return false;/* point_segment_idx -- contains -1 */
    else
        return true;/* point_segment_idx -- does not contain -1 */
}
void  MaSegProcess::remaining_idx_remove_idx(long long int id) {
    std::vector<long long int>::iterator it = std::find((*remaining_idx).begin(), (*remaining_idx).end(), id);
    if (it == (*remaining_idx).end())
        std::cout << "remaining_idx_remove_idx -- error " << std::endl;
    int index = std::distance((*remaining_idx).begin(), it);
    (*remaining_idx).erase((*remaining_idx).begin() + index);
}

inline size_t MaSegProcess::findseed8r(float seed_radius_thres,floatList *ma_radius) {
    for (size_t i = 1; i < this->size; ++i) {
        if (point_segment_idx[i] != -1 && (*ma_radius)[i] < seed_radius_thres)
            return i;
    }
}

inline size_t MaSegProcess::findseed(size_t initial_seed_idx) {
    for (size_t i = initial_seed_idx; i < this->size; ++i) {
        if (point_segment_idx[i] == -1)
            return i;
    }
}
size_t MaSegProcess::findseed_random() {
    std::random_device rd; // obtain a random number from hardware
    std::mt19937 eng(rd()); // seed the generator
    std::uniform_int_distribution<> distr(0, (*remaining_idx).size()-1); // define the range
    auto random_id = distr(eng);
    return (*remaining_idx)[random_id];
}
inline bool MaSegProcess::valid_candidate_bisec(float bisec_thres,size_t idx1, size_t idx2, ma_Geometry &maGeometry) {
    Vector bisec1 = maGeometry.ma_bisector[idx1];
    Vector bisec2 = maGeometry.ma_bisector[idx2];
    /*
    std::cout << "idx1--"<< idx1<<"  idx2--"<< idx2<<"  dot product : "<<bisec1 * bisec2 << '\n';
    std::cout << "bisec_thres -- " << bisec_thres << "\n";
    std::cout << (bisec1 * bisec2 > bisec_thres) << "\n";
    */
    return bisec1 * bisec2 > bisec_thres;
}

inline bool MaSegProcess::valid_candidate_spokecross(float spokecross_thres, size_t idx1, size_t idx2, ma_Geometry &maGeometry) {
    Vector crosnorm1 = maGeometry.ma_normal[idx1];
    Vector crosnorm2 = maGeometry.ma_normal[idx2];
    return crosnorm1 * crosnorm2 > spokecross_thres;
}

inline bool MaSegProcess::valid_candidate_balloverlap(float balloverlap_thres, size_t idx1, size_t  idx2, mat_data &madata) {
    Point p1 = madata.ma_coords[idx1];
    Point p2 = madata.ma_coords[idx2];
    auto dis = Vrui::Geometry::dist(p1, p2);
    auto r1 = madata.ma_radius[idx1];
    auto r2 = madata.ma_radius[idx2];
    return (r1 + r2) / dis > balloverlap_thres;

}

bool MaSegProcess::validateCandidate(MaSeg_power &power,size_t idx1, size_t idx2,
    mat_data &madata, ma_Geometry &maGeometry) {
    if (power.method == masb::bisector)
        return valid_candidate_bisec(power.bisec_thres,idx1,idx2, maGeometry);
    if (power.method == masb::spokecross)
        return valid_candidate_spokecross(power.spokecross_thres, idx1, idx2, maGeometry);
    if (power.method == masb::balloverlap)
        return valid_candidate_balloverlap(power.balloverlap_thres, idx1, idx2, madata);
    if (power.method == masb::combinBisecAndSpcros)
        return valid_candidate_bisec(power.bisec_thres, idx1, idx2, maGeometry)&&
        valid_candidate_spokecross(power.spokecross_thres, idx1, idx2, maGeometry);
    if (power.method == masb::combinBallAndSpcros)
        return valid_candidate_balloverlap(power.balloverlap_thres, idx1, idx2, madata)&&
        valid_candidate_spokecross(power.spokecross_thres, idx1, idx2, maGeometry);

    else
        std::cout << "not having this segmentation method yet, plz check.\n";
}

void MaSegProcess::grow_sheet(MaSeg_power &power,size_t initial_seed_idx,
    mat_data &madata, ma_Geometry &maGeometry) {
    stack<size_t> seeds;
    size_tList idx_in_sheet;
    int in_out_flag = madata.in_out[initial_seed_idx];
    /*
    std::cout << "initial_seed_idx--" << initial_seed_idx << "\n";
    std::cout << "bisector-thresh--" << power.bisec_thres << "\n";
    std::cout << "in_out_flag -- " << in_out_flag << "\n";
    */
    seeds.push(initial_seed_idx);
    idx_in_sheet.push_back(initial_seed_idx);
    point_segment_idx[initial_seed_idx] = sheet_counter;
    //remaining_idx_remove_idx(initial_seed_idx);

    while (!seeds.empty()) {
        size_t seed_idx = seeds.top();
        seeds.pop();
        Point seed_pt = madata.ma_coords[seed_idx];
        kdtree2::KDTreeResultVector neighbours;
        madata.kdtree_ma_coords->n_nearest(seed_pt, power.k_neib, neighbours);
        for (auto &candidate : neighbours) {
            /*
            std::cout << "\ncandidate.idx -- " << candidate.idx << '\t';
            std::cout << "point_segment_idx -- " << point_segment_idx[candidate.idx] << "\t";
            std::cout << "candidate in_out -- " << madata.in_out[candidate.idx] << "\n";
            */
            if (point_segment_idx[candidate.idx] == -1 &&
                madata.in_out[candidate.idx] == in_out_flag &&
                validateCandidate(power,seed_idx, candidate.idx, madata, maGeometry)
                ) {
                //std::cout << "candidate validate, push-- " << candidate.idx << " into seeds,sheet,segment\n";
                seeds.push(candidate.idx);
                idx_in_sheet.push_back(candidate.idx);
                point_segment_idx[candidate.idx] = sheet_counter;
                //std::cout << "point_segment_idx[candidate.idx] -- " << point_segment_idx[candidate.idx] << "\n";
            }
        }
    }
    if (idx_in_sheet.size() < power.mincount) {
        for (size_t &idx : idx_in_sheet) {
            point_segment_idx[idx] = -1;
        }
        point_segment_idx[initial_seed_idx] = 0;
    }
    else {
        shape.push_back(idx_in_sheet);
        shape_inout.push_back(in_out_flag);
        sheet_counter++;
        //for (size_t &idx : idx_in_sheet) {
        //    remaining_idx_remove_idx(idx);
        //}
    }
}

void MaSegProcess::processing(MaSeg_power &power,mat_data &madata,ma_Geometry &maGeometry) {
    //power.update();
    std::cout << "MaSegProcess::processing -- min num is " << power.mincount << std::endl;
    this->size = madata.ma_ptsize;
    point_segment_idx.resize(size, -1);

    //std::vector<long long int> v(size); // vector with 100 ints.
    //std::iota(std::begin(v), std::end(v), 0); // Fill with 0, 1, ..., 99.
    //*remaining_idx = v;

    if (madata.kdtree_ma_coords == NULL)
        madata.kdtree_ma_coords = new kdtree2::KDTree(madata.ma_coords, true);
    madata.kdtree_ma_coords->sort_results = true;
    size_t initial_seed_idx = 0;
    //size_t initial_seed_idx = findseed_random();
    while (1) {
        grow_sheet(power,initial_seed_idx,madata, maGeometry);
        if (if_all_segmented())
            break;
        else
            initial_seed_idx = findseed(initial_seed_idx);
            //initial_seed_idx = findseed_random();
    }
    std::cout << "segmentation finished\n";
}

