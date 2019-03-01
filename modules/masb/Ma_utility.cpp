#include "Ma_utility.h"
#include <iostream>

using namespace masb;
void idx_filter::processing(ma_data &madata, ma_Geometry &maGeometry,
    intList &remaining_idx, mat_data &remainingData, intList &remainingma_in_out, 
    ma_Geometry &remainingGeometry) {

    std::vector<int> ma_is_interior(madata.m * 2, 0);
    std::fill_n(ma_is_interior.begin(), madata.m, 1);
        
    remainingData.m = madata.m;
    remainingData.coords = madata.coords;
    remainingData.normals = madata.normals;
    size_t size = remaining_idx.size();
    PointList re_coords_(size);
    floatList re_radius_(size);
    intList re_qidx_(size);
    VectorList re_bisector_(size);
    floatList re_SeperationAng_(size);
    remainingma_in_out.reserve(size);

    for (int &i : remaining_idx) {
        re_coords_.push_back((*madata.ma_coords)[i]);
        re_radius_.push_back((*madata.ma_radius)[i]);
        re_qidx_.push_back(madata.ma_qidx[i]);//why different????
        re_bisector_.push_back((*maGeometry.ma_bisector)[i]);
        re_SeperationAng_.push_back((*maGeometry.ma_SeperationAng)[i]);
        remainingma_in_out.push_back(ma_is_interior[i]);
    }
    remainingData.ma_coords = re_coords_;
    remainingData.ma_radius = re_radius_;
    remainingData.ma_qidx = re_qidx_;
    remainingGeometry.ma_bisector = &re_bisector_;
    remainingGeometry.ma_SeperationAng = &re_SeperationAng_;
    std::cout << "the number of remaining ma_point -- " << remaining_idx.size() << "\n";
}

void idx_filter::processing(ma_data &madata,
    intList &remaining_idx, mat_data &remainingData, intList &remainingma_in_out) {

    std::vector<int> ma_is_interior(madata.m * 2, 0);
    std::fill_n(ma_is_interior.begin(), madata.m, 1);

    remainingData.m = remaining_idx.size();
    PointList re_coords_(remainingData.m);
    floatList re_radius_(remainingData.m);
    intList re_qidx_(remainingData.m);
    VectorList re_bisector_(remainingData.m);
    floatList re_SeperationAng_(remainingData.m);
    remainingma_in_out.reserve(remainingData.m);

    for (int i : remaining_idx) {
        re_coords_.push_back((*madata.ma_coords)[i]);
        re_radius_.push_back((*madata.ma_radius)[i]);
        re_qidx_.push_back(madata.ma_qidx[i]);//why different????
        //re_bisector_.push_back((*maGeometry.ma_bisector)[i]);
        //re_SeperationAng_.push_back((*maGeometry.ma_SeperationAng)[i]);
        remainingma_in_out.push_back(ma_is_interior[i]);
    }
    remainingData.ma_coords = re_coords_;
    remainingData.ma_radius = re_radius_;
    remainingData.ma_qidx = re_qidx_;
    //remainingGeometry.ma_bisector = &re_bisector_;
    //remainingGeometry.ma_SeperationAng = &re_SeperationAng_;
}

void MaPt_in_oneTrace::processing(PointList &ma_coords, ma_Geometry &maGeometry, pt_Traces &traces) {

}