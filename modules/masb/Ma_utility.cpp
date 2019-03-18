#include "Ma_utility.h"
#include <iostream>

using namespace masb;
void idx_filter::processing(ma_data &madata, ma_Geometry &maGeometry,
    intList &remaining_idx, mat_data &remainingData,
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
    intList re_in_out_; re_in_out_.reserve(size);
    //remainingma_in_out.reserve(size);

    for (auto &idx : remaining_idx) {
        re_coords_.push_back((*madata.ma_coords)[idx]);
        re_radius_.push_back((*madata.ma_radius)[idx]);
        re_qidx_.push_back(madata.ma_qidx[idx]);
        re_bisector_.push_back(maGeometry.ma_bisector[idx]);
        re_SeperationAng_.push_back(maGeometry.ma_SeperationAng[idx]);
        re_in_out_.push_back(ma_is_interior[idx]);
    }
    remainingData.ma_coords = re_coords_;
    remainingData.ma_radius = re_radius_;
    remainingData.ma_qidx = re_qidx_;
    remainingData.in_out = re_in_out_;
    remainingGeometry.ma_bisector = re_bisector_;
    remainingGeometry.ma_SeperationAng = re_SeperationAng_;

    auto it = std::find(re_in_out_.begin(), re_in_out_.end(), 0);
    remainingData.in_ptsize = std::distance(re_in_out_.begin(), it);
    remainingData.out_ptsize = remainingData.ma_ptsize - remainingData.in_ptsize;
    std::cout << "the number of remaining ma_point -- " << remaining_idx.size() << "\n";
}
void idx_filter::processing(mat_data &madata, ma_Geometry &maGeometry,
    size_tList &remaining_idx, mat_data &remainingData,
    ma_Geometry &remainingGeometry) {
    remainingData.m = madata.m;

    size_t size = remaining_idx.size();
    remainingData.ma_ptsize = size;
    PointList re_coords_;    re_coords_.reserve(size);
    floatList re_radius_;    re_radius_.reserve(size);
    intList re_qidx_;    re_qidx_.reserve(size);
    VectorList re_bisector_;    re_bisector_.reserve(size);
    floatList re_SeperationAng_;    re_SeperationAng_.reserve(size);
    intList re_in_out_; re_in_out_.reserve(size);

    for (auto &idx : remaining_idx) {
        re_coords_.push_back(madata.ma_coords[idx]);
        re_radius_.push_back(madata.ma_radius[idx]);
        re_qidx_.push_back(madata.ma_qidx[idx]);
        re_in_out_.push_back(madata.in_out[idx]);
        re_bisector_.push_back(maGeometry.ma_bisector[idx]);
        re_SeperationAng_.push_back(maGeometry.ma_SeperationAng[idx]);
    }
    remainingData.ma_coords = re_coords_;
    remainingData.ma_radius = re_radius_;
    remainingData.ma_qidx = re_qidx_;
    remainingData.in_out = re_in_out_;
    remainingGeometry.ma_bisector = re_bisector_;
    remainingGeometry.ma_SeperationAng = re_SeperationAng_;
        
    remainingData.in_ptsize = std::count(remainingData.in_out.begin(), 
        remainingData.in_out.end(), 1);
    remainingData.out_ptsize = remainingData.ma_ptsize - remainingData.in_ptsize;
}
void idx_filter::processing(ma_data &madata, intList &remaining_idx,
    mat_data &remainingData) {

    
}

void splitInOut::processing(mat_data &madata, ma_Geometry &maGeometry, bool InFlag,
    mat_data &madata_new, ma_Geometry &maGeometry_new) {
    if (InFlag == true) {
        size = madata.in_ptsize;
        begin = 0;
        madata_new.in_ptsize = size;
        madata_new.out_ptsize = 0;
        madata_new.in_out.resize(size, 1);
    }
    else {
        size = madata.out_ptsize;
        begin = madata.in_ptsize;
        madata_new.in_ptsize = 0;
        madata_new.out_ptsize = size;
        madata_new.in_out.resize(size, 0);
    }
    madata_new.ma_ptsize = size;
    madata_new.ma_coords = VectorSlice(madata.ma_coords, size, begin);
    madata_new.ma_radius = VectorSlice(madata.ma_radius, size, begin);
    madata_new.ma_qidx = VectorSlice(madata.ma_qidx, size, begin);
    maGeometry_new.ma_bisector = VectorSlice(maGeometry.ma_bisector, size, begin);
    maGeometry_new.ma_SeperationAng = VectorSlice(maGeometry.ma_SeperationAng,size, begin);
}

