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
inline size_t MaSegProcess::findseed(size_t initial_seed_idx) {
    for (size_t i = initial_seed_idx; i < this->size; ++i) {
        if (point_segment_idx[i] == -1)
            return i;
    }
}
/*
void  MaSegProcess::remaining_idx_remove_idx(long long int id) {
    std::vector<long long int>::iterator it = std::find((*remaining_idx).begin(), (*remaining_idx).end(), id);
    if (it == (*remaining_idx).end())
        std::cout << "remaining_idx_remove_idx -- error " << std::endl;
    int index = std::distance((*remaining_idx).begin(), it);
    (*remaining_idx).erase((*remaining_idx).begin() + index);
}
size_t MaSegProcess::findseed_random() {
    std::random_device rd; // obtain a random number from hardware
    std::mt19937 eng(rd()); // seed the generator
    std::uniform_int_distribution<> distr(0, (*remaining_idx).size()-1); // define the range
    auto random_id = distr(eng);
    return (*remaining_idx)[random_id];
}
*/
inline bool MaSegProcess::valid_candidate_bisec(float bisec_thres,size_t idx1, size_t idx2, MAT&mat) {
    Vector bisec1 = mat.bisector[idx1];
    Vector bisec2 = mat.bisector[idx2];
    /*
    std::cout << "idx1--"<< idx1<<"  idx2--"<< idx2<<"  dot product : "<<bisec1 * bisec2 << '\n';
    std::cout << "bisec_thres -- " << bisec_thres << "\n";
    std::cout << (bisec1 * bisec2 > bisec_thres) << "\n";
    */
    return bisec1 * bisec2 > bisec_thres;
}
inline bool MaSegProcess::valid_candidate_sepAng(float sepAng_thres, size_t idx1, size_t  idx2, MAT&mat) {
    auto diff = mat.seperationAng[idx1] - mat.seperationAng[idx2];
    if (diff < 0)
        diff = -diff;
    return diff < sepAng_thres;
}
inline bool MaSegProcess::valid_candidate_spokecross(float spokecross_thres, size_t idx1, size_t idx2, MAT&mat) {
    Vector crosnorm1 = mat.ma_direction[idx1];
    Vector crosnorm2 = mat.ma_direction[idx2];
    return crosnorm1 * crosnorm2 > spokecross_thres;
}

inline bool MaSegProcess::valid_candidate_balloverlap(float balloverlap_thres, size_t idx1, size_t  idx2, MAT&mat) {
    Point p1 = mat.atom[idx1];
    Point p2 = mat.atom[idx2];
    auto dis = Vrui::Geometry::dist(p1, p2);
    auto r1 = mat.radius[idx1];
    auto r2 = mat.radius[idx2];
    return (r1 + r2) / dis > balloverlap_thres;

}

bool MaSegProcess::validateCandidate(MaSeg_power &power,size_t idx1, size_t idx2, MAT&mat) {
    if (power.method == masb::bisector)
        return valid_candidate_bisec(power.bisec_thres,idx1,idx2, mat);
    if (power.method == masb::spokecross)
        return valid_candidate_spokecross(power.spokecross_thres, idx1, idx2, mat);
    if (power.method == masb::balloverlap)
        return valid_candidate_balloverlap(power.balloverlap_thres, idx1, idx2, mat);
    if (power.method == masb::combinBisecAndSpcros)
        return valid_candidate_bisec(power.bisec_thres, idx1, idx2, mat)&&
        valid_candidate_spokecross(power.spokecross_thres, idx1, idx2, mat);
    if (power.method == masb::combinBallAndSpcros)
        return valid_candidate_balloverlap(power.balloverlap_thres, idx1, idx2, mat)&&
        valid_candidate_spokecross(power.spokecross_thres, idx1, idx2, mat);

    else
        std::cout << "not having this segmentation method yet, plz check.\n";
}

void MaSegProcess::grow_sheet(MaSeg_power &power, size_t initial_seed_idx, 
    MAT &mat, kdtree2::KDTree *kdtree_ma_atoms) {
    stack<size_t> seeds;
    size_tList idx_in_sheet;
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
        Point seed_pt = mat.atom[seed_idx];
        kdtree2::KDTreeResultVector neighbours;
        kdtree_ma_atoms->n_nearest(seed_pt, power.k_neib, neighbours);
        for (auto &candidate : neighbours) {
            /*
            std::cout << "\ncandidate.idx -- " << candidate.idx << '\t';
            std::cout << "point_segment_idx -- " << point_segment_idx[candidate.idx] << "\t";
            std::cout << "candidate in_out -- " << madata.in_out[candidate.idx] << "\n";
            */
            if (point_segment_idx[candidate.idx] == -1 &&
                validateCandidate(power,seed_idx, candidate.idx, mat)) {
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
        sheet_counter++;
        //for (size_t &idx : idx_in_sheet) {
        //    remaining_idx_remove_idx(idx);
        //}
    }
}

void MaSegProcess::processing(MaSeg_power &power, MAT &mat){
    std::cout << "MaSegProcess::processing -- min num is " << power.mincount << std::endl;
    this->size = mat.matSize;
    point_segment_idx.resize(size, -1);

    //std::vector<long long int> v(size); // vector with 100 ints.
    //std::iota(std::begin(v), std::end(v), 0); // Fill with 0, 1, ..., 99.
    //*remaining_idx = v;

    auto kdtree_ma_atoms = new kdtree2::KDTree(mat.atom, true);
    kdtree_ma_atoms->sort_results = true;
    size_t initial_seed_idx = 0;
    //size_t initial_seed_idx = findseed_random();
    while (1) {
        grow_sheet(power,initial_seed_idx,mat, kdtree_ma_atoms);
        if (if_all_segmented())
            break;
        else
            initial_seed_idx = findseed(initial_seed_idx);
            //initial_seed_idx = findseed_random();
    }
    std::cout << "segmentation finished\n";
}

