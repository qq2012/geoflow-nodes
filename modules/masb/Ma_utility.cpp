#include "Ma_utility.h"
#include <iostream>

using namespace masb;
void idx_filter::processing(ma_data &madata, ma_Geometry &maGeometry,
    intList &remaining_idx, mat_data &remainingData, intList &remainingma_in_out, 
    ma_Geometry &remainingGeometry) {

    std::vector<int> ma_is_interior(madata.m * 2, 0);
    std::fill_n(ma_is_interior.begin(), madata.m, 1);

    remainingData.m = madata.m;

    size_t size = remaining_idx.size();
    remainingData.ma_ptsize = size;
    PointList re_coords_;    re_coords_.reserve(size);
    floatList re_radius_;    re_radius_.reserve(size);
    intList re_qidx_;    re_qidx_.reserve(size);
    VectorList re_bisector_;    re_bisector_.reserve(size);
    floatList re_SeperationAng_;    re_SeperationAng_.reserve(size);

    remainingma_in_out.reserve(size);

    for (int &idx : remaining_idx) {
        re_coords_.push_back((*madata.ma_coords)[idx]);
        re_radius_.push_back((*madata.ma_radius)[idx]);
        re_qidx_.push_back(madata.ma_qidx[idx]);
        re_bisector_.push_back(maGeometry.ma_bisector[idx]);
        re_SeperationAng_.push_back(maGeometry.ma_SeperationAng[idx]);
        remainingma_in_out.push_back(ma_is_interior[idx]);
    }
    remainingData.ma_coords = re_coords_;
    remainingData.ma_radius = re_radius_;
    remainingData.ma_qidx = re_qidx_;
    remainingGeometry.ma_bisector = re_bisector_;
    remainingGeometry.ma_SeperationAng = re_SeperationAng_;
    std::cout << "the number of remaining ma_point -- " << remaining_idx.size() << "\n";
}

void idx_filter::processing(ma_data &madata,
    intList &remaining_idx, mat_data &remainingData, intList &remainingma_in_out) {

    
}

void MaPt_in_oneTrace::processing(PointList &ma_coords, ma_Geometry &maGeometry, pt_Traces &traces) {

}