#include "Ma_Segmentation_processing.h"
#include <iostream>
#include <cmath>
#include <stack>
using namespace masb;
using namespace std;

//static const double PI = 3.14159265358979323846264338327950288;

void MaSeg_power::update() {
    this->bisec_thres = cos((this->bisec_thres / 180.0)*PI);
    this->bisecavg_thres = cos((this->bisecavg_thres / 180.0)*PI);
    this->bisecdiff_thres = cos((this->bisecdiff_thres / 180.0)*PI);
    this->theta_thres = this->theta_thres / 180.0 *PI;
    this->spokecross_thres = cos((this->spokecross_thres / 180.0)*PI);
    /*
    self.p_bisecthres = math.cos((self.p['bisec_thres'] / 180.0) * math.pi)
    self.p_bisecavgthres = math.cos((self.p['bisecavg_thres'] / 180.0) * math.pi)
    self.p_bisecdiffthres = math.cos((self.p['bisecdiff_thres'] / 180.0) * math.pi)
    self.p_normalthres = math.cos((5.0 / 180.0) * math.pi)
    self.p_thetathres_1 = (self.p['theta_thres'] / 180.0) * math.pi # during bisect growing
    self.p_thetathres_2 = (self.p['theta_thres'] / 180.0) * math.pi # during theta growing
    self.p_k = self.p['k']
    self.p_balloverlap_thres = self.p['balloverlap_thres']
    self.p_mincount = self.p['mincount']
    self.p_spokecross_thres = math.cos((self.p['spokecross_thres'] / 180.0) * math.pi)
   */
}

inline bool MaSegProcess::if_all_segmented() {
    if (std::find(point_segment_idx.begin(), point_segment_idx.end(), -1) != point_segment_idx.end())
        return false;/* point_segment_idx -- contains -1 */
    else
        return true;/* point_segment_idx -- does not contain -1 */
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

bool MaSegProcess::validateCandidate(MaSeg_power &power,size_t idx1, size_t idx2,
    mat_data &madata, ma_Geometry &maGeometry) {
    if (power.method == masb::bisector)
        return valid_candidate_bisec(power.bisec_thres,idx1,idx2, maGeometry);
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
    }
}

void MaSegProcess::processing(MaSeg_power &power,mat_data &madata,ma_Geometry &maGeometry) {
    power.update();
    this->size = madata.ma_ptsize;
    point_segment_idx.resize(size, -1);
    if (madata.kdtree_ma_coords == NULL)
        madata.kdtree_ma_coords = new kdtree2::KDTree(madata.ma_coords, true);
    madata.kdtree_ma_coords->sort_results = true;
    size_t initial_seed_idx = 0;
    while (1) {
        grow_sheet(power,initial_seed_idx,madata, maGeometry);
        if (if_all_segmented())
            break;
        else
            initial_seed_idx = findseed(initial_seed_idx);
    }
    std::cout << "segmentation finished\n";
}

