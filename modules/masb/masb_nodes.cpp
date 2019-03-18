#include "masb_nodes.hpp"
#include <iostream>
#include <algorithm>
#include <string>

namespace geoflow::nodes::mat {

void ComputeMedialAxisNode::process(){
  auto point_collection = input("points").get<PointCollection>();
  auto normals_vec3f = input("normals").get<vec3f>();

  masb::ma_parameters params;
  params.initial_radius = param<float>("initial_radius");
  params.denoise_preserve = param<double>("denoise_preserve");
  params.denoise_planar = param<double>("denoise_planar");
  params.nan_for_initr = param<bool>("nan_for_initr");

  masb::ma_data madata;
  madata.m = point_collection.size();
  
  masb::PointList coords;
  coords.reserve(madata.m);
  for(auto& p : point_collection) {
    coords.push_back(masb::Point(p.data()));
  }
  masb::VectorList normals;
  normals.reserve(madata.m);
  for(auto& n : normals_vec3f) {
    normals.push_back(masb::Vector(n.data()));
  }
  masb::PointList ma_coords_(madata.m*2);
  std::vector<float> ma_radius_(madata.m*2);
  std::vector<int> ma_qidx_(madata.m*2);
  
  
  madata.coords = &coords;
  madata.normals = &normals;
  madata.ma_coords = &ma_coords_;
  madata.ma_radius = &ma_radius_;
  madata.ma_qidx = ma_qidx_.data();

  masb::compute_masb_points(params, madata);

  vec1i ma_qidx;
  ma_qidx.reserve(madata.m*2);
  for(size_t i=0 ; i<madata.m*2; ++i) {
    ma_qidx.push_back(madata.ma_qidx[i]);
  }

  PointCollection ma_coords;
  ma_coords.reserve(madata.m*2);
  for(auto& c : *madata.ma_coords) {
    ma_coords.push_back({c[0], c[1], c[2]});
  }

  vec1f ma_radius;
  ma_radius.reserve(madata.m * 2);
  for (auto& r : *madata.ma_radius) {
      ma_radius.push_back(r);
  }

  vec1i ma_is_interior(madata.m*2, 0);
  std::fill_n(ma_is_interior.begin(), madata.m, 1);

  output("coords_masb").set(coords);
  output("ma_coords").set(ma_coords);
  output("ma_qidx").set(ma_qidx);
  output("ma_radius").set(ma_radius);
  output("ma_is_interior").set(ma_is_interior);
}


void ComputeNormalsNode::process(){
  auto point_collection = input("points").get<PointCollection>();

  masb::normals_parameters params;
  params.k = param<int>("K");

  masb::ma_data madata;
  madata.m = point_collection.size();
  masb::PointList coords;
  coords.reserve(madata.m);
  for(auto& p : point_collection) {
    coords.push_back(masb::Point(p.data()));
  }
  masb::VectorList normals(madata.m);
  madata.coords = &coords;
  madata.normals = &normals;

  masb::compute_normals(params, madata);

  vec3f normals_vec3f;
  normals_vec3f.reserve(madata.m);
  for(auto& n : *madata.normals) {
    normals_vec3f.push_back({n[0], n[1], n[2]});
  }

  output("normals").set(normals_vec3f);
}

void MaGeometryNode::process() {
    auto point_collection = input("points").get<PointCollection>();
    auto normals_vec3f = input("normals").get<vec3f>();
    auto ma_coords_collection = input("ma_coords").get<PointCollection>();
    auto ma_qidx_vec1i = input("ma_qidx").get<vec1i>();

    masb::ma_data madata;
    madata.m = point_collection.size();

    masb::PointList coords;
    coords.reserve(madata.m);
    for (auto& p : point_collection) {
        coords.push_back(masb::Point(p.data()));
    }
    masb::VectorList normals;
    normals.reserve(madata.m);
    for (auto& n : normals_vec3f) {
        normals.push_back(masb::Vector(n.data()));
    }    
    masb::PointList ma_coords;
    ma_coords.reserve(madata.m * 2);
    for (auto& p : ma_coords_collection) {
        ma_coords.push_back(masb::Point(p.data()));
    }
    masb::intList ma_qidx;
    ma_qidx.reserve(madata.m * 2);
    for (size_t i = 0; i < madata.m * 2; ++i) {
        ma_qidx.push_back(ma_qidx_vec1i[i]);
    }
    vec1i ma_is_interior(madata.m * 2, 0);
    std::fill_n(ma_is_interior.begin(), madata.m, 1);
    madata.coords = &coords;
    madata.normals = &normals;
    madata.ma_coords = &ma_coords;
    madata.ma_qidx = &(ma_qidx[0]);//?????????????????????????????????????

    std::vector<float> ma_SeparationAng_(madata.m * 2);
    masb::VectorList ma_bisector_(madata.m * 2);
    masb::ma_Geometry maGeometry;
    maGeometry.ma_SeperationAng = ma_SeparationAng_;
    maGeometry.ma_bisector = ma_bisector_;

    masb::compute_ma_geometry(madata, maGeometry);

    vec1f ma_SeparationAng;
    ma_SeparationAng.reserve(madata.m * 2);
    for (auto& c : maGeometry.ma_SeperationAng) {
        ma_SeparationAng.push_back(c);
    }
    vec3f ma_bisector;
    ma_bisector.reserve(madata.m * 2);
    for (auto& n : maGeometry.ma_bisector) {
        ma_bisector.push_back({ n[0], n[1], n[2] });
    }

    output("ma_SeparationAng").set(ma_SeparationAng);
    output("ma_bisector").set(ma_bisector);
}

void FilterRNode::process() {
    auto ma_radius_vec1f = input("ma_radius").get<vec1f>();

    float radius = param<float>("filterRadius");

    vec1i remaining_idx;
    for (size_t i = 0; i< ma_radius_vec1f.size()-1; ++i) {
        if (ma_radius_vec1f[i] < radius)
            remaining_idx.push_back(i);
    }
    std::cout << remaining_idx.data() << '\n';
    output("remaining_idx").set(remaining_idx);
}

void MedialSegmentNode::process() {
    auto remaining_idx_vec1i = input("remaining_idx").get<vec1i>();
    auto ma_coords_collection = input("ma_coords").get<PointCollection>();
    auto ma_qidx_vec1i = input("ma_qidx").get<vec1i>();
    auto ma_radius_vec1f = input("ma_radius").get<vec1f>();
    auto ma_SeparationAng_vec1f = input("ma_SeparationAng").get<vec1f>();
    auto ma_bisector_vec3f = input("ma_bisector").get<vec3f>();

    masb::MaSeg_power params;
    params.mincount = param<int>("mincount");
    params.maxcount = param<int>("maxcount");
    params.method = this->current_method;

    switch (this->current_method) {
    case masb::bisector:
        params.bisec_thres = param<float>("bisec_thres");
        std::cout << "method -- bisector\n"; 
        break;
    case masb::radius:
        params.balloverlap_thres = param<float>("balloverlap_thres"); break;
    case masb::thirdopt:
        std::cout << "not ready\n"; break;
    }

    
    masb::intList remaining_idx;
    for (auto& i : remaining_idx_vec1i) {
        remaining_idx.push_back(int(i));
    }

    masb::ma_data madata;
    masb::mat_data remainingData;
    madata.m = ma_coords_collection.size() / 2;
    masb::PointList ma_coords;
    ma_coords.reserve(madata.m*2);
    for (auto& p : ma_coords_collection) {
        ma_coords.push_back(masb::Point(p.data()));
    }
    masb::floatList ma_radius;
    ma_radius.reserve(madata.m * 2);
    for (auto& r : ma_radius_vec1f) {
        ma_radius.push_back(r);
    }
    masb::intList ma_qidx;
    ma_qidx.reserve(madata.m * 2);
    for (int & i : ma_qidx_vec1i) {
        ma_qidx.push_back(i);
    }
    madata.ma_coords = &ma_coords;
    madata.ma_radius = &ma_radius;
    madata.ma_qidx = &(ma_qidx[0]);//?????????????????????????
    
    masb::ma_Geometry maGeometry, remainingGeometry;
    masb::floatList ma_SeparationAng;
    ma_SeparationAng.reserve(madata.m * 2);
    for (auto & a : ma_SeparationAng_vec1f) {
        ma_SeparationAng.push_back(float(a));
    }
    masb::VectorList ma_bisector;
    ma_bisector.reserve(madata.m * 2);
    for (auto& v : ma_bisector_vec3f) {
        ma_bisector.push_back(masb::Vector(v.data()));
    }
    maGeometry.ma_SeperationAng = ma_SeparationAng;
    maGeometry.ma_bisector = ma_bisector;

    //masb::intList remainingma_in_out;
    masb::idx_filter filter;
    filter.processing(madata, maGeometry, remaining_idx, 
        remainingData, remainingGeometry);

    masb::mat_data madata_in, madata_out;
    masb::ma_Geometry maGeometry_in, maGeometry_out;
    masb::splitInOut splier;
    splier.processing(remainingData, remainingGeometry, true, madata_in, maGeometry_in);
    splier.processing(remainingData, remainingGeometry, false, madata_out, maGeometry_out);

    masb::MaSegProcess segmentation;
    segmentation.processing(params, madata_in, maGeometry_in);
    auto seg_in_ = segmentation.point_segment_idx;
    auto seg_all_ = segmentation.point_segment_idx;
    auto sheet_in_ = segmentation.shape;
    auto sheet_all_ = segmentation.shape;

    segmentation.processing(params, madata_out, maGeometry_out);
    auto seg_out_ = segmentation.point_segment_idx;
    auto sheet_out_ = segmentation.shape;

    seg_all_.insert(seg_all_.end(), segmentation.point_segment_idx.begin(), segmentation.point_segment_idx.end());
    sheet_all_.insert(sheet_all_.end(), segmentation.shape.begin(), segmentation.shape.end());

    vec1i seg_idx_;
    for (auto & idx : seg_all_)
        seg_idx_.push_back(idx);

    PointCollection ma_coords_;
    ma_coords_.reserve(remainingData.ma_ptsize);
    for (auto& c : remainingData.ma_coords) {
        ma_coords_.push_back({ c[0], c[1], c[2] });
    }

    output("seg_id").set(seg_idx_);
    output("sheet_all").set(sheet_all_);
    output("ma_coords").set(ma_coords_);
    output("madata_in").set(madata_in);
    output("madata_out").set(madata_out);
    output("maGeometry_in").set(maGeometry_in);
    output("maGeometry_out").set(maGeometry_out);
    output("seg_in").set(seg_in_);
    output("seg_out").set(seg_out_);
    output("sheet_in").set(sheet_in_);
    output("sheet_out").set(sheet_out_);

}

void MaPt_in_oneTraceNode::process() {
    auto madata = input("ma_coords").get<masb::mat_data>();
    auto maGeometry = input("maGeometry").get<masb::ma_Geometry>();
    auto sheets = input("sheets").get<masb::Sheet_idx_List>();

    masb::pt_Trace_pram  params;
    params.SearchRadius = param<float>("SearchRadius");
    masb::MaPt_in_oneTrace tracer;
    tracer.processing(params,madata, maGeometry, sheets);

}

void ExtractCandidatePtNode::process() {


}
void ConnectCandidatePtNode::process() {

}
}

