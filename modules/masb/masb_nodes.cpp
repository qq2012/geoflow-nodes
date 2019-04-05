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

  vec1f ma_radii(madata.m*2);
  for(size_t i=0; i<madata.m*2; ++i) {
    ma_radii.push_back( Vrui::Geometry::dist(coords[i%madata.m], ma_coords_[i]) );
  }

  vec1i ma_is_interior(madata.m*2, 0);
  std::fill_n(ma_is_interior.begin(), madata.m, 1);

  output("coords_masb").set(coords);
  output("ma_coords").set(ma_coords);
  output("ma_radii").set(ma_radii);
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
    
    masb::floatList ma_SeparationAng_(madata.m * 2);
    masb::VectorList ma_bisector_(madata.m * 2);
    masb::VectorList ma_normal_(madata.m * 2);
    
    masb::ma_Geometry maGeometry;
    maGeometry.ma_SeperationAng = ma_SeparationAng_;
    maGeometry.ma_bisector = ma_bisector_;
    maGeometry.ma_normal = ma_normal_;
    
    masb::compute_ma_geometry(madata, maGeometry);

    vec1f ma_SeparationAng;
    ma_SeparationAng.reserve(madata.m * 2);
    for (auto& c : maGeometry.ma_SeperationAng) {
        ma_SeparationAng.push_back(c);
    }
    vec3f ma_bisector;
    ma_bisector.reserve(madata.m * 2);
    for (auto n : maGeometry.ma_bisector) {
        ma_bisector.push_back({ n[0], n[1], n[2] });
    }
    vec3f ma_normal;
    ma_normal.reserve(madata.m * 2);
    for (auto& n : maGeometry.ma_normal) {
        ma_normal.push_back({ n[0], n[1], n[2] });
    }

    LineStringCollection bisec;
    bisec.reserve(madata.ma_coords->size());
    size_t i = 0;
    for (auto &p : *madata.ma_coords) {
        LineString temp;
        LineString cp, cq;
        if (madata.ma_qidx[i] != -1) {
            temp.push_back({ p[0],p[1],p[2] });
            auto end = p + maGeometry.ma_bisector[i];
            temp.push_back({ end[0],end[1],end[2] });
            /*
            cp.push_back({ p[0],p[1],p[2] });
            cq.push_back({ p[0],p[1],p[2] });
            if (i >= madata.m)
                auto p = (*madata.coords)[i - madata.m];
            else
                auto p = (*madata.coords)[i];
            auto q = (*madata.coords)[madata.ma_qidx[i]];
            cp.push_back({ p[0],p[1],p[2] });
            cq.push_back({ q[0],q[1],q[2] });
            */
        }
        bisec.push_back(temp);
        //bisec.push_back(cp);
        //bisec.push_back(cq);
        i++;
    }

    output("ma_SeparationAng").set(ma_SeparationAng);
    output("ma_bisector").set(ma_bisector);
    output("ma_normal").set(ma_normal);
    output("bisec").set(bisec);
}

void FilterRNode::process() {
    auto ma_radius_vec1f = input("ma_radius").get<vec1f>();

    float radius = param<float>("filterRadius");

    vec1i remaining_idx;
    for (size_t i = 0; i< ma_radius_vec1f.size()-1; ++i) {
        if (ma_radius_vec1f[i] < radius)
            remaining_idx.push_back(i);
    }
    //std::cout << remaining_idx.data() << '\n';
    output("remaining_idx").set(remaining_idx);
}

void MedialSegmentNode::process() {
    auto remaining_idx_vec1i = input("remaining_idx").get<vec1i>();
    auto ma_coords_collection = input("ma_coords").get<PointCollection>();
    auto ma_qidx_vec1i = input("ma_qidx").get<vec1i>();
    auto ma_radius_vec1f = input("ma_radius").get<vec1f>();
    auto ma_SeparationAng_vec1f = input("ma_SeparationAng").get<vec1f>();
    auto ma_bisector_vec3f = input("ma_bisector").get<vec3f>();
    auto ma_normal_vec3f = input("ma_normal").get<vec3f>();
    
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
    masb::VectorList ma_normal;
    ma_normal.reserve(madata.m * 2);
    for (auto& v : ma_normal_vec3f) {
        ma_normal.push_back(masb::Vector(v.data()));
    }
    maGeometry.ma_SeperationAng = ma_SeparationAng;
    maGeometry.ma_bisector = ma_bisector;
    maGeometry.ma_normal = ma_normal;

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
    auto madata = input("madata").get<masb::mat_data>();
    auto maGeometry = input("maGeometry").get<masb::ma_Geometry>();
    auto sheets = input("sheets").get<masb::Sheet_idx_List>();

    masb::pt_Trace_pram  params;
    params.SearchRadius = param<float>("SearchRadius");
    masb::MaPt_in_oneTrace tracer;
    tracer.processing(params,madata, maGeometry, sheets);

    vec1i seg_idx_;
    PointCollection candidate_r_, candidate_cos_;
    candidate_r_.reserve(tracer.candidate_size);
    candidate_cos_.reserve(tracer.candidate_size);

    int i = 0;
    for (auto canInSheet : tracer.candidate_r) {
        for (auto p : canInSheet) {
            candidate_r_.push_back({ p[0], p[1], p[2] });
        }
        vec1i temp;
        temp.resize(canInSheet.size(), i++);//check should start from 1;
        seg_idx_.insert(seg_idx_.end(), temp.begin(), temp.end());
    }
    for (auto canInSheet : tracer.candidate_cos) {
        for (auto p : canInSheet) {
            candidate_cos_.push_back({ p[0], p[1], p[2] });
        }
    }
    LineStringCollection traces;
    //int j =0;
    for (auto &traces_in_1_sheet : tracer.all_traces) {
        for (auto & a_trace : traces_in_1_sheet) {
            LineString a_trace_LineString;
            a_trace_LineString.reserve(a_trace.size());
            for (auto p : a_trace) {
                a_trace_LineString.push_back({ p[0], p[1], p[2] });
            }
            traces.push_back(a_trace_LineString);
        }
    }

    output("candidate_r").set(candidate_r_);
    output("candidate_cos").set(candidate_cos_);
    output("lable").set(seg_idx_);
    output("traces").set(traces);// TT_line_string_collection);
    
}

void ExtractCandidatePtNode::process() {
    auto madata = input("madata").get<masb::mat_data>();
    auto maGeometry = input("maGeometry").get<masb::ma_Geometry>();
    auto sheets = input("sheets").get<masb::Sheet_idx_List>();

    masb::ExtractCandidatePt_pram  params;
    params.SearchRadius = param<float>("SearchRadius");
    masb::ExtractCandidatePt extractor;
    extractor.processing(params, madata, maGeometry, sheets);

    PointCollection candidate_r_, candidate_cos_;
    candidate_r_.reserve(extractor.candidate_size);
    candidate_cos_.reserve(extractor.candidate_size);

    int i = 0;
    for (auto canInSheet : extractor.candidate_r) {
        for (auto p : canInSheet) {
            candidate_r_.push_back({ p[0], p[1], p[2] });
        }
    }
    for (auto canInSheet : extractor.candidate_cos) {
        for (auto p : canInSheet) {
            candidate_cos_.push_back({ p[0], p[1], p[2] });
        }
    }
    output("candidate_r").set(candidate_r_);
    output("candidate_cos").set(candidate_cos_);
}
void ReadCandidatePtNode::process() {
    PointCollection candidate_r_;
    vec3f direction_;
    vec1i seg_id_;

    //std::vector<double> v;
    std::ifstream infile;
    std::string filepath;
    filepath = (std::string) "C:/Users/wangq/Downloads/thesis/p3_data/2cpp.ply";
    infile.open(filepath);
    std::string dummyLine;
    int size;
    int i = 0;
    while (i < 12) {
        std::getline(infile, dummyLine);
        if (i == 3) {
            auto elems = masb::split2nums(dummyLine, ' ');
            size = ::atof(elems[2].c_str());
        }
        i++;
    }
    candidate_r_.reserve(size);
    direction_.reserve(size);
    seg_id_.reserve(size);

    std::string numbers;
    while (!infile.eof()) {
        std::getline(infile, numbers);
        auto elems = masb::split2nums(numbers, ' ');
        float x = ::atof(elems[0].c_str());
        float y = ::atof(elems[1].c_str());
        float z = ::atof(elems[2].c_str());
        int seg_idx = ::atof(elems[3].c_str());
        float nx = ::atof(elems[4].c_str());
        float ny = ::atof(elems[5].c_str());
        float nz = ::atof(elems[6].c_str());
        candidate_r_.push_back({ x,y,z });
        seg_id_.push_back(seg_idx);
        direction_.push_back({ nx,ny,nz });
        //i++;
    }
    infile.close();
    if (candidate_r_.size() != size)
        std::cout << "candidate error" << std::endl;
    if (seg_id_.size() != size)
        std::cout << "candidate error" << std::endl;
    if (direction_.size() != size)
        std::cout << "candidate error" << std::endl;

    
    output("candidate_r").set(candidate_r_);
    output("directon").set(direction_);
    output("seg_id").set(seg_id_);
}

void ReadCandidatePtWithBisecNode::process() {
    PointCollection candidate_r_;
    vec3f direction_, bisector_p_, bisector_q_;
    vec1i seg_id_;

    std::ifstream infile;
    std::string filepath;
    //filepath = (std::string) "C:/Users/wangq/Downloads/thesis/p3_data/candidatept_r45SepAng15_notAlongBisec_midz_all_WithBisector.ply";
    filepath = (std::string) "C:/Users/wangq/Downloads/thesis/p3_data/candidatept_r45SepAng10_alongBisec_sheet3_8_WithBisector.ply";
    infile.open(filepath);
    std::string dummyLine;
    int size;
    int i = 0;
    while (i < 18) {
        std::getline(infile, dummyLine);
        if (i == 3) {
            auto elems = masb::split2nums(dummyLine, ' ');
            size = ::atof(elems[2].c_str());
        }
        i++;
    }
    candidate_r_.reserve(size);
    direction_.reserve(size);
    seg_id_.reserve(size);
    bisector_p_.reserve(size);
    bisector_q_.reserve(size);

    std::string numbers;
    while (!infile.eof()) {
        std::getline(infile, numbers);
        auto elems = masb::split2nums(numbers, ' ');
        float x = ::atof(elems[0].c_str());
        float y = ::atof(elems[1].c_str());
        float z = ::atof(elems[2].c_str());
        int seg_idx = ::atof(elems[3].c_str());
        float nx = ::atof(elems[4].c_str());
        float ny = ::atof(elems[5].c_str());
        float nz = ::atof(elems[6].c_str());
        float bpx = ::atof(elems[7].c_str());
        float bpy = ::atof(elems[8].c_str());
        float bpz = ::atof(elems[9].c_str());
        float bqx = ::atof(elems[10].c_str());
        float bqy = ::atof(elems[11].c_str());
        float bqz = ::atof(elems[12].c_str());
        candidate_r_.push_back({ x,y,z });
        seg_id_.push_back(seg_idx);
        direction_.push_back({ nx,ny,nz });
        bisector_p_.push_back({ bpx,bpy,bpz });
        bisector_q_.push_back({ bqx,bqy,bqz });
        //i++;
    }
    infile.close();
    if (candidate_r_.size() != size)
        std::cout << "candidate error" << std::endl;
    if (seg_id_.size() != size)
        std::cout << "candidate error" << std::endl;
    if (direction_.size() != size)
        std::cout << "candidate error" << std::endl;
    if (bisector_p_.size() != size)
        std::cout << "bisector_p error" << std::endl;
    if (bisector_q_.size() != size)
        std::cout << "bisector_q_ error" << std::endl;

    output("candidate_r").set(candidate_r_);
    output("directon").set(direction_);
    output("seg_id").set(seg_id_);
    output("bisector_p").set(bisector_p_);
    output("bisector_q").set(bisector_q_);
}

void ConnectCandidatePtNode::process() {

    auto pointCloud_ptcollection = input("pointCloud").get<PointCollection>();
    auto candidate_r_ptcollection = input("candidate").get<PointCollection>();
    auto directon_vec3f = input("directon").get<vec3f>();
    auto seg_id_vec1i = input("seg_id").get<vec1i>();
    auto bisector_p_vec3f = input("bisector_p").get<vec3f>(); 
    auto bisector_q_vec3f = input("bisector_q").get<vec3f>();

    masb::PointList pointCloud;
    pointCloud.reserve(pointCloud_ptcollection.size());
    for (auto& p : pointCloud_ptcollection) {
        pointCloud.push_back(masb::Point(p.data()));
    }

    size_t size = candidate_r_ptcollection.size();
    masb::PointList candidate_r;
    candidate_r.reserve(size);
    for (auto& p : candidate_r_ptcollection) {
        candidate_r.push_back(masb::Point(p.data()));
    }
    masb::VectorList directon;
    directon.reserve(size);
    for (auto& n : directon_vec3f) {
        auto n_normlize = masb::Vector(n.data()).normalize();
        directon.push_back(n_normlize);
    }
    masb::VectorList bisector_p;
    bisector_p.reserve(size);
    for (auto& bp : bisector_p_vec3f) {
        auto bp_normlize = masb::Vector(bp.data()).normalize();
        bisector_p.push_back(bp_normlize);
    }
    masb::VectorList bisector_q;
    bisector_q.reserve(size);
    for (auto& bq : bisector_q_vec3f) {
        auto bq_normlize = masb::Vector(bq.data()).normalize();
        bisector_q.push_back(bq_normlize);
    }
    masb::intList seg_id;
    seg_id.reserve(size);
    for (auto id : seg_id_vec1i) {
        seg_id.push_back(id);
    }
    masb::intList filter, filter2;
    ridge::segment segmentList;
    masb::intList idList, idList2;
    ridge::line symple_segmentList;
    ridge::line line_segmentList;
    ridge::line smoothLine;
    masb::intList symple_idList;
    
    //ridge::connectCandidatePt8Spline(pointCloud, candidate_r, seg_id,
    //    filter2, line_segmentList, idList2);

    ridge::connectCandidatePt8MST(pointCloud, candidate_r, seg_id,
        filter, segmentList, idList, symple_segmentList, symple_idList);

    ridge::connectCandidatePtSmooth(symple_segmentList, smoothLine);

    vec1i filter_;
    filter_.reserve(filter.size());
    for (auto i : filter)
        filter_.push_back(i);

    LineStringCollection segment_vis_;
    segment_vis_.reserve(segmentList.size());
    for (auto &s : segmentList) {
        LineString tmp;
        tmp.push_back({ s.first[0],s.first[1],s.first[2] });
        tmp.push_back({ s.second[0],s.second[1] ,s.second[2] });
        segment_vis_.push_back(tmp);
    }
    vec1i idList_;
    idList_.reserve(idList.size());
    for (auto i : idList)
        idList_.push_back(i);

    LineStringCollection longest_seg_vis_;
    longest_seg_vis_.reserve(symple_segmentList.size());
    for (auto &a_line : symple_segmentList) {
        LineString tmp;
        for (auto &p : a_line) {
            tmp.push_back({ p[0],p[1],p[2] });
        }
        longest_seg_vis_.push_back(tmp);
    }
    vec1i symple_idList_;
    symple_idList_.reserve(symple_idList.size());
    for (auto i : symple_idList)
        symple_idList_.push_back(i);

    LineStringCollection smoothLine_vis_;
    smoothLine_vis_.reserve(smoothLine.size());
    for (auto& a_smooth_line : smoothLine) {
        LineString tmp;
        for (auto &p : a_smooth_line) {
            tmp.push_back({ p[0],p[1],p[2] });
        }
        smoothLine_vis_.push_back(tmp);
    }

    LineStringCollection directon_vis_;
    LineStringCollection directon2_vis_;
    LineStringCollection bisec_p_vis_;
    LineStringCollection bisec_q_vis_;
    for (int j = 0; j < size; ++j) {
        auto p = candidate_r[j];
        auto v = directon[j];
        auto bp = bisector_p[j];
        auto bq = bisector_q[j];
        LineString tmp, tmp2, tmp_bp,tmp_bq;
        tmp.push_back({ p[0],p[1],p[2] });
        tmp2.push_back({ p[0],p[1],p[2] });
        tmp_bp.push_back({ p[0],p[1],p[2] });
        tmp_bq.push_back({ p[0],p[1],p[2] });
        auto end = p + v;
        auto v2= Vrui::Geometry::cross(bp, bq);
        v2 = v2.normalize();
        auto end2 = p + v2;
        auto end_bp = p + bp;
        auto end_bq = p + bq;
        tmp.push_back({ end[0] ,end[1] ,end[2] });
        tmp2.push_back({ end2[0] ,end2[1] ,end2[2] });
        tmp_bp.push_back({ end_bp[0] ,end_bp[1] ,end_bp[2] });
        tmp_bq.push_back({ end_bq[0] ,end_bq[1] ,end_bq[2] });
        directon_vis_.push_back(tmp);
        directon2_vis_.push_back(tmp2);
        bisec_p_vis_.push_back(tmp_bp);
        bisec_q_vis_.push_back(tmp_bq);
    }
    output("filter").set(filter_);
    output("ridge").set(segment_vis_);
    output("ridgeId").set(idList_);
    output("directon_vis").set(directon_vis_);//directon_vis
    output("directon2_vis").set(directon2_vis_);
    output("bisector_p_vis").set(bisec_p_vis_);
    output("bisector_q_vis").set(bisec_q_vis_);
    output("longest_path").set(longest_seg_vis_);
    output("longest_id").set(symple_idList_);
    output("smoothLine").set(smoothLine_vis_);

    
}
void PLYLoaderNode::process() {
    int k = param<int>("thinning_factor");
    std::string filepath = param<std::string>("filepath");
    std::cout << "start load point cloud\n";
    std::ifstream infile;
    infile.open(filepath);
    std::string dummyLine;
    int size;
    int i = 0;
    while (i < 8) {
        std::getline(infile, dummyLine);
        //std::cout << i << "--" << dummyLine << std::endl;
        if (i == 3) {
            auto elems = masb::split2nums(dummyLine, ' ');
            size = ::atof(elems[2].c_str());
        }
        i++;
    }
    PointCollection PointCloud_;
    PointCloud_.reserve(size / k + 1);
    std::string xyz;
    for (int j = 0; j < size ; j++) {
        std::getline(infile, xyz);
        if (j%k == 0) {
            auto elems = masb::split2nums(xyz, ' ');
            float x = ::atof(elems[0].c_str());
            float y = ::atof(elems[1].c_str());
            float z = ::atof(elems[2].c_str());
            PointCloud_.push_back({ x,y,z });
        }
    }
    output("PointCloud").set(PointCloud_);
}
}

