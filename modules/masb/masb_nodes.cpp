#include "masb_nodes.hpp"
#include <iostream>
#include <algorithm>
#include <string>

namespace geoflow::nodes::mat {

    void ComputeMedialAxisNode::process() {
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
        for (auto& p : point_collection) {
            coords.push_back(masb::Point(p.data()));
        }
        masb::VectorList normals;
        normals.reserve(madata.m);
        for (auto& n : normals_vec3f) {
            normals.push_back(masb::Vector(n.data()));
        }
        masb::PointList ma_coords_(madata.m * 2);
        std::vector<float> ma_radius_(madata.m * 2);
        std::vector<int> ma_qidx_(madata.m * 2);


        madata.coords = &coords;
        madata.normals = &normals;
        madata.ma_coords = &ma_coords_;
        madata.ma_radius = &ma_radius_;
        madata.ma_qidx = ma_qidx_.data();

        masb::compute_masb_points(params, madata);

        vec1i ma_qidx;
        ma_qidx.reserve(madata.m * 2);
        for (size_t i = 0; i < madata.m * 2; ++i) {
            ma_qidx.push_back(madata.ma_qidx[i]);
        }

        PointCollection ma_coords;
        ma_coords.reserve(madata.m * 2);
        for (auto& c : *madata.ma_coords) {
            ma_coords.push_back({ c[0], c[1], c[2] });
        }

        vec1f ma_radius;
        ma_radius.reserve(madata.m * 2);
        for (auto& r : *madata.ma_radius) {
            ma_radius.push_back(r);
        }

        vec1f ma_radii(madata.m * 2);
        for (size_t i = 0; i < madata.m * 2; ++i) {
            ma_radii.push_back(Vrui::Geometry::dist(coords[i%madata.m], ma_coords_[i]));
        }

        vec1i ma_is_interior(madata.m * 2, 0);
        std::fill_n(ma_is_interior.begin(), madata.m, 1);

        //output("coords_masb").set(coords);
        output("ma_coords").set(ma_coords);
        output("ma_radii").set(ma_radii);
        output("ma_radius").set(ma_radius);
        output("ma_qidx").set(ma_qidx);
        output("ma_is_interior").set(ma_is_interior);
    }

    void ComputeNormalsNode::process() {
        auto point_collection = input("points").get<PointCollection>();

        masb::normals_parameters params;
        params.k = param<int>("K");

        masb::ma_data madata;
        madata.m = point_collection.size();
        masb::PointList coords;
        coords.reserve(madata.m);
        for (auto& p : point_collection) {
            coords.push_back(masb::Point(p.data()));
        }
        masb::VectorList normals(madata.m);
        madata.coords = &coords;
        madata.normals = &normals;

        masb::compute_normals(params, madata);

        vec3f normals_vec3f;
        normals_vec3f.reserve(madata.m);
        for (auto& n : *madata.normals) {
            normals_vec3f.push_back({ n[0], n[1], n[2] });
        }

        output("normals").set(normals_vec3f);
    }

    void MaGeometryNode::process() {
        auto point_cloud = input("points").get<PointCollection>();
        auto normals_vec3f = input("normals").get<vec3f>();
        auto ma_coords_collection = input("ma_coords").get<PointCollection>();
        auto ma_radius_vec1f = input("ma_radius").get<vec1f>();
        auto ma_qidx_vec1i = input("ma_qidx").get<vec1i>();
        auto isInterior = param<bool>("isInterior");

        masb::ma_data madata;
        madata.m = point_cloud.size();

        masb::PointList coords;
        coords.reserve(madata.m);
        for (auto& p : point_cloud) {
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
        masb::floatList ma_radius;
        ma_radius.reserve(madata.m * 2);
        for (auto r : ma_radius_vec1f) {
            ma_radius.push_back(r);
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
        madata.ma_radius = &ma_radius;
        madata.ma_qidx = &(ma_qidx[0]);//?????????????????????????????????????

        masb::MAT mat;
        mat.isInterior = isInterior;

        masb::compute_ma_geometry_new(madata, mat, isInterior);

        mat.matSize = mat.atom.size();
        PointCollection atom_;atom_.reserve(mat.matSize);
        for (auto &p : mat.atom)
            atom_.push_back({ p[0],p[1],p[2] });
        //vec1f radius_; radius_.reserve(mat.matSize);
        //for (auto r : mat.radius)
        //    radius_.push_back(r);
        PointCollection unshrikingGroundPoint_;unshrikingGroundPoint_.reserve(madata.m - mat.matSize);
        for (auto &p : mat.unshrikingGroundPoint)
            unshrikingGroundPoint_.push_back({ p[0],p[1],p[2] });
        //mat.fp = fp_;//feature point p
        //mat.fq = fq_;
        vec3f sp_reverse_norm_;sp_reverse_norm_.reserve(mat.matSize);
        for (auto &v : mat.sp_reverse_norm)
            sp_reverse_norm_.push_back({ v[0],v[1],v[2] });
        vec3f sp_;sp_.reserve(mat.matSize);//spoken vector
        for (auto &v : mat.sp)
            sp_.push_back({ v[0],v[1],v[2] });
        //mat.sq = sq_;
        vec3f bisector_; bisector_.reserve(mat.matSize);
        for (auto &v : mat.bisector)
            bisector_.push_back({ v[0],v[1],v[2] });
        //mat.seperationAng = seperationAng_;
        vec3f ma_direction_; ma_direction_.reserve(mat.matSize);
        for (auto &v : mat.ma_direction)
            ma_direction_.push_back({ v[0],v[1],v[2] });

        output("atom").set(atom_);
        output("unshrikingGroundPoint").set(unshrikingGroundPoint_);
        output("sp_reverse_norm").set(sp_reverse_norm_);
        output("spokeVectorP").set(sp_);
        //add_output("spokeVectorQ", typeid(vec3f));
        output("bisector").set(bisector_);
        output("ma_direction").set(ma_direction_);
        output("mat").set(mat);
    }

    void FilterRNode::process() {
        auto ma_radius_vec1f = input("ma_radius").get<vec1f>();

        float radius = param<float>("filterRadius");

        vec1i remaining_idx;
        for (size_t i = 0; i < ma_radius_vec1f.size() - 1; ++i) {
            if (ma_radius_vec1f[i] < radius)
                remaining_idx.push_back(i);
        }
        //std::cout << remaining_idx.data() << '\n';
        output("remaining_idx").set(remaining_idx);
    }

    void MedialSegmentNode::process() {
        auto mat = input("mat").get<masb::MAT>();
        
        masb::MaSeg_power params;
        params.mincount = param<int>("mincount");
        params.maxcount = param<int>("maxcount");
        params.method = this->current_method;

        //static const double PI = 3.14159265358979323846264338327950288;
        switch (this->current_method) {
        case masb::bisector:
            params.bisec_thres = cos((param<float>("bisec_thres") / 180.0)*PI);
            std::cout << "method -- bisector\n"
                << "bisec_thres " << param<float>("bisec_thres") << "degree\n"
                << "cos(bisec_thres) " << params.bisec_thres << std::endl;
            break;
        case masb::spokecross:
            params.spokecross_thres = cos((param<float>("spokecross_thres") / 180.0)*PI);
            std::cout << "method -- spokecross\n"
                << "spokecross_thres " << param<float>("spokecross_thres") << "degree\n"
                << "cos(spokecross_thres) " << params.spokecross_thres << std::endl;
            break;
        case masb::balloverlap:
            params.balloverlap_thres = param<float>("balloverlap_thres");
            std::cout << "method -- balloverlap\n"
                << "balloverlap__thres " << params.balloverlap_thres << " scalar."
                << std::endl;
            break;
        case masb::combinBisecAndSpcros:
            params.bisec_thres = cos((param<float>("bisec_thres") / 180.0)*PI);
            params.spokecross_thres = cos((param<float>("spokecross_thres") / 180.0)*PI);
            std::cout << "method -- combine bisector and spokecross\n"
                << "bisec_thres " << param<float>("bisec_thres") << "degree\n"
                << "cos(bisec_thres) " << params.bisec_thres
                << "spokecross_thres " << param<float>("spokecross_thres") << "degree\n"
                << "cos(spokecross_thres) " << params.spokecross_thres
                << std::endl;
            break;
        case masb::combinBallAndSpcros:
            params.balloverlap_thres = param<float>("balloverlap_thres");
            params.spokecross_thres = cos((param<float>("spokecross_thres") / 180.0)*PI);
            std::cout << "method -- combine balloverlap and spokecross\n"
                << "balloverlap__thres " << params.balloverlap_thres << " scalar\n"
                << "spokecross_thres " << param<float>("spokecross_thres") << "degree\n"
                << "cos(spokecross_thres) " << params.spokecross_thres
                << std::endl;
            break;
        }
        
        masb::MaSegProcess segmentation;
        segmentation.processing(params, mat);

        vec1i seg_id_;
        seg_id_.reserve(mat.matSize);
        for (auto& i : segmentation.point_segment_idx) {
            seg_id_.push_back((int)i);
        }

        output("seg_id").set(seg_id_);

    }

    void MaPt_in_oneTraceNode::process() {
        auto madata = input("madata").get<masb::mat_data>();
        auto maGeometry = input("maGeometry").get<masb::ma_Geometry>();
        auto sheets = input("sheets").get<masb::Sheet_idx_List>();

        masb::pt_Trace_pram  params;
        params.SearchRadius = param<float>("SearchRadius");
        masb::MaPt_in_oneTrace tracer;
        tracer.processing(params, madata, maGeometry, sheets);

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

        auto pointcloud = input("point cloud coords").get<PointCollection>();
        auto pointcloud_norm = input("point cloud normal").get<vec3f>();

        masb::ExtractCandidatePt_pram  params;
        params.SearchRadius = param<float>("SearchRadius");
        params.deviationAng_thres = cos((param<float>("deviationAng_thres") / 180.0)*PI);
        std::cout << "deviationAng_thres is " << param<float>("deviationAng_thres") << " degree"
            << "cos (deviationAng_thres) is" << params.deviationAng_thres << std::endl;

        masb::PointCloud PointCloud;
        masb::PointList pc_coords_;
        pc_coords_.reserve(madata.m);
        for (auto &p : pointcloud)
            pc_coords_.push_back(masb::Point(p.data()));
        masb::VectorList pc_normals_;
        pc_normals_.reserve(madata.m);
        for (auto &v : pointcloud_norm)
            pc_normals_.push_back(masb::Vector(v.data()));
        PointCloud.coords = pc_coords_;
        PointCloud.normals = pc_normals_;


        masb::ExtractCandidatePt extractor;
        extractor.processing(params, madata, maGeometry, sheets, PointCloud);

        PointCollection candidate_r_, candidate_cos_;
        PointCollection candidate_bisector_avg_;
        vec1i seg_id_;
        vec3f spokV_cp_, spokV_cq_, direction_;

        int s = extractor.can_pt_r.size();
        candidate_r_.reserve(s);
        candidate_cos_.reserve(s);
        candidate_bisector_avg_.reserve(s);
        seg_id_.reserve(s);
        spokV_cp_.reserve(s);
        spokV_cq_.reserve(s);
        direction_.reserve(s);
        for (int i = 0; i < s; ++i) {
            if (extractor.seg_id[i] == 0)
                continue;
            auto c1 = extractor.can_pt_r[i];
            candidate_r_.push_back({ c1[0],c1[1],c1[2] });

            auto c2 = extractor.can_pt_cos[i];
            candidate_cos_.push_back({ c2[0],c2[1,c2[2]] });

            auto c3 = extractor.can_pt_bisector_avg[i];
            candidate_bisector_avg_.push_back({ c3[0],c3[1],c3[2] });

            seg_id_.push_back(extractor.seg_id[i]);

            auto v1 = extractor.cp[i];
            spokV_cp_.push_back({ v1[0],v1[1],v1[2] });

            auto v2 = extractor.cq[i];
            spokV_cq_.push_back({ v2[0],v2[1],v2[2] });

            auto d = extractor.direction[i];
            direction_.push_back({ d[0],d[1],d[2] });
        }

        std::cout << "ExtractCandidatePtNode:: there are " << candidate_r_.size() << " candidate points" << std::endl;

        output("candidate_r").set(candidate_r_);
        output("candidate_cos").set(candidate_cos_);
        output("candidate_bisector_avg").set(candidate_bisector_avg_);
        output("seg_id").set(seg_id_);
        output("spoke_cp").set(spokV_cp_);
        output("spoke_cq").set(spokV_cq_);
        output("direction").set(direction_);

        /*
        PointCollection candidate_r_, candidate_cos_;
        candidate_r_.reserve(extractor.candidate_size);
        candidate_cos_.reserve(extractor.candidate_size);
        int i = 0;
        for (auto canInSheet : extractor.candidate_r) {
            if (i == 0) {
                i++;
                continue;
            }
            for (auto p : canInSheet) {
                candidate_r_.push_back({ p[0], p[1], p[2] });
            }
            i++;
        }
        int j = 0;
        for (auto canInSheet : extractor.candidate_cos) {
            if (j == 0) {
                j++;
                continue;
            }
            for (auto p : canInSheet) {
                candidate_cos_.push_back({ p[0], p[1], p[2] });
            }
            j++;
        }
        output("candidate_r").set(candidate_r_);
        output("candidate_cos").set(candidate_cos_);
        */
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

        std::string filepath = param<std::string>("filepath");
        std::cout << "start load point cloud\n";

        PointCollection candidate_r_;
        vec3f direction_, bisector_p_, bisector_q_;
        vec1i seg_id_;

        std::ifstream infile;
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
            float cpx = ::atof(elems[7].c_str());
            float cpy = ::atof(elems[8].c_str());
            float cpz = ::atof(elems[9].c_str());
            float cqx = ::atof(elems[10].c_str());
            float cqy = ::atof(elems[11].c_str());
            float cqz = ::atof(elems[12].c_str());

            masb::Vector cp = { cpx,cpy,cpz };
            cp.normalize();
            masb::Vector cq = { cqx,cqy,cqz };
            cq.normalize();
            auto b = cp + cq;
            b.normalize();
            auto n = Vrui::Geometry::cross(cp, cq);
            n.normalize();

            candidate_r_.push_back({ x,y,z });
            seg_id_.push_back(seg_idx);
            direction_.push_back({ n[0],n[1],n[2] });
            bisector_p_.push_back({ cp[0],cp[1],cp[2] });
            bisector_q_.push_back({ cp[0],cp[1],cp[2] });
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
    void ReadSegmentRestltNode::process() {
        std::string filepath = param<std::string>("filepath");
        std::cout << "start load Segment Restlt\n";

        PointCollection ma_coords_;
        vec1i seg_id_;
        vec1f radius_, theta_;
        vec3f cp_, cq_, bisector_, direction_;
        LineStringCollection direction_vis_;

        std::ifstream infile;
        infile.open(filepath);
        std::string dummyLine;
        int size;
        int i = 0;
        while (i < 17) {
            std::getline(infile, dummyLine);
            if (i == 3) {
                auto elems = masb::split2nums(dummyLine, ' ');
                size = ::atof(elems[2].c_str());
            }
            i++;
        }

        ma_coords_.reserve(size);
        seg_id_.reserve(size);
        radius_.reserve(size);
        theta_.reserve(size);
        cp_.reserve(size);
        cq_.reserve(size);
        bisector_.reserve(size);
        direction_.reserve(size);
        direction_vis_.reserve(size);

        std::string numbers;
        while (!infile.eof()) {
            std::getline(infile, numbers);
            auto elems = masb::split2nums(numbers, ' ');
            float x = ::atof(elems[0].c_str());
            float y = ::atof(elems[1].c_str());
            float z = ::atof(elems[2].c_str());
            int seg_id = ::atof(elems[3].c_str());
            float r = ::atof(elems[4].c_str());
            float cpx = ::atof(elems[5].c_str());
            float cpy = ::atof(elems[6].c_str());
            float cpz = ::atof(elems[7].c_str());
            float cqx = ::atof(elems[8].c_str());
            float cqy = ::atof(elems[9].c_str());
            float cqz = ::atof(elems[10].c_str());
            float theta = ::atof(elems[11].c_str());



            masb::Vector cp = { cpx,cpy,cpz };
            cp.normalize();
            masb::Vector cq = { cqx,cqy,cqz };
            cq.normalize();
            auto b = cp + cq;
            b.normalize();
            auto n = Vrui::Geometry::cross(cp, cq);
            n.normalize();

            ma_coords_.push_back({ x,y,z });
            seg_id_.push_back(seg_id);
            radius_.push_back(r);
            theta_.push_back(theta);
            cp_.push_back({ cp[0],cp[1],cp[2] });
            cq_.push_back({ cq[0],cq[1],cq[2] });
            bisector_.push_back({ b[0],b[1],b[2] });
            direction_.push_back({ n[0],n[1],n[2] });
            LineString temp;
            temp.push_back({ x,y,z });
            temp.push_back({ x + n[0],y + n[1],z + n[2] });
            direction_vis_.push_back({ temp });
            //i++;
        }
        infile.close();
        if (ma_coords_.size() != size)
            std::cout << "ma_coords error" << std::endl;
        if (seg_id_.size() != size)
            std::cout << "candidate error" << std::endl;
        if (radius_.size() != size)
            std::cout << "radius error" << std::endl;
        if (theta_.size() != size)
            std::cout << "radius error" << std::endl;
        if (direction_.size() != size)
            std::cout << "candidate error" << std::endl;
        if (bisector_.size() != size)
            std::cout << "bisector error" << std::endl;
        if (cp_.size() != size)
            std::cout << "bisector_p error" << std::endl;
        if (cq_.size() != size)
            std::cout << "bisector_q_ error" << std::endl;

        output("ma_coords").set(ma_coords_);
        output("seg_id").set(seg_id_);
        output("radius").set(radius_);
        output("cp").set(cp_);
        output("cq").set(cq_);
        output("bisector").set(bisector_);
        output("directon").set(direction_);
        output("theta").set(theta_);
        output("direction_vis").set(direction_vis_);
    }
    void ReadJunctionPtNode::process() {

        std::string filepath = param<std::string>("filepath");
        std::cout << "start load JunctionPt\n";

        PointCollection JunctionPt_;
        ridge::int_pair_vec sheet_sheet_;

        std::ifstream infile;
        infile.open(filepath);
        std::string dummyLine;
        int size;
        int i = 0;
        while (i < 10) {
            std::getline(infile, dummyLine);
            if (i == 3) {
                auto elems = masb::split2nums(dummyLine, ' ');
                size = ::atof(elems[2].c_str());
            }
            i++;
        }
        JunctionPt_.reserve(size);
        sheet_sheet_.reserve(size);

        std::string numbers;
        while (!infile.eof()) {
            std::getline(infile, numbers);
            auto elems = masb::split2nums(numbers, ' ');
            float x = ::atof(elems[0].c_str());
            float y = ::atof(elems[1].c_str());
            float z = ::atof(elems[2].c_str());
            int i = ::atof(elems[3].c_str());
            int j = ::atof(elems[4].c_str());

            JunctionPt_.push_back({ x,y,z });
            sheet_sheet_.push_back(std::make_pair(i, j));
            //int_pair(i,j)

            //i++;
        }
        infile.close();
        if (JunctionPt_.size() != size)
            std::cout << "Junction Point error" << std::endl;
        if (sheet_sheet_.size() != size)
            std::cout << "sheet_sheet error" << std::endl;

        output("ma_coords").set(JunctionPt_);
        output("sheet-sheet").set(sheet_sheet_);
    }
    void ReadAdjacencyNode::process() {
        std::string filepath = param<std::string>("filepath");
        std::cout << "start load Adjacency map\n";

        ridge::int_pair_vec sheet_sheet_;

        std::ifstream infile;
        infile.open(filepath);
        std::string dummyLine;
        int size;
        int i = 0;
        while (i < 7) {
            std::getline(infile, dummyLine);
            if (i == 3) {
                auto elems = masb::split2nums(dummyLine, ' ');
                size = ::atof(elems[2].c_str());
            }
            i++;
        }
        sheet_sheet_.reserve(size);

        std::string numbers;
        while (!infile.eof()) {
            std::getline(infile, numbers);
            auto elems = masb::split2nums(numbers, ' ');
            int i = ::atof(elems[0].c_str());
            int j = ::atof(elems[1].c_str());
            sheet_sheet_.push_back(std::make_pair(i, j));
            //int_pair(i,j)
            //i++;
        }
        infile.close();
        if (sheet_sheet_.size() != size)
            std::cout << "sheet_sheet error" << std::endl;

        output("sheet-sheet").set(sheet_sheet_);
    }
    void ConnectCandidatePtNode::process() {

        auto pointCloud_ptcollection = input("pointCloud").get<PointCollection>();
        auto candidate_r_ptcollection = input("candidate").get<PointCollection>();
        auto directon_vec3f = input("directon").get<vec3f>();
        auto seg_id_vec1i = input("seg_id").get<vec1i>();
        auto bisector_p_vec3f = input("bisector_p").get<vec3f>();
        auto bisector_q_vec3f = input("bisector_q").get<vec3f>();
        auto adjacency = input("adjacency").get<ridge::int_pair_vec>();

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
        ridge::segment segmentList, segmentList_nosegid;
        masb::intList idList, idList2;
        ridge::line symple_segmentList;
        ridge::line line_segmentList;
        ridge::line smoothLine;
        masb::intList symple_idList;

        //ridge::connectCandidatePt8Spline(pointCloud, candidate_r, seg_id,
        //    filter2, line_segmentList, idList2);

        std::cout << "STARTING METHOD 1 -- with SEG\n";
        ridge::connectCandidatePt8MST(pointCloud, candidate_r, seg_id,
            filter, segmentList, idList, symple_segmentList, symple_idList);

        //ridge::connectCandidatePtSmooth(symple_segmentList, smoothLine);

        ridge::line smoothLines = symple_segmentList;
        ridge::FindTopology(smoothLines, symple_idList, adjacency);

        //std::cout << "STARTING METHOD 2 -- NO SEG-ID\n";
        //ridge::connectCandidatePt8MST_nosegid(pointCloud, candidate_r, filter2, segmentList_nosegid);

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

        LineStringCollection segment_nosegid_vis_;
        segment_nosegid_vis_.reserve(segmentList_nosegid.size());
        for (auto &s : segmentList_nosegid) {
            LineString tmp;
            tmp.push_back({ s.first[0],s.first[1],s.first[2] });
            tmp.push_back({ s.second[0],s.second[1] ,s.second[2] });
            segment_nosegid_vis_.push_back(tmp);
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
        /*
        PointCollection smoothpoint_vis_;
        LineStringCollection smoothLine_vis_;
        smoothLine_vis_.reserve(smoothLine.size());
        for (auto& a_smooth_line : smoothLine) {
            LineString tmp;
            for (auto &p : a_smooth_line) {
                tmp.push_back({ p[0],p[1],p[2] });
                smoothpoint_vis_.push_back({ p[0],p[1],p[2] });
            }
            smoothLine_vis_.push_back(tmp);
        }
        */
        LineStringCollection topology_vis_;
        topology_vis_.reserve(smoothLines.size());
        for (auto& a_line : smoothLines) {
            LineString tmp;
            for (auto &p : a_line) {
                tmp.push_back({ p[0],p[1],p[2] });
            }
            topology_vis_.push_back(tmp);
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
            LineString tmp, tmp2, tmp_bp, tmp_bq;
            tmp.push_back({ p[0],p[1],p[2] });
            tmp2.push_back({ p[0],p[1],p[2] });
            tmp_bp.push_back({ p[0],p[1],p[2] });
            tmp_bq.push_back({ p[0],p[1],p[2] });
            auto end = p + v;
            auto v2 = Vrui::Geometry::cross(bp, bq);
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
        output("ridge_nosegid_vis_").set(segment_nosegid_vis_);
        output("ridgeId").set(idList_);
        output("directon_vis").set(directon_vis_);//directon_vis
        output("directon2_vis").set(directon2_vis_);
        output("bisector_p_vis").set(bisec_p_vis_);
        output("bisector_q_vis").set(bisec_q_vis_);
        output("longest_path").set(longest_seg_vis_);
        output("longest_id").set(symple_idList_);
        //output("smoothLine").set(smoothLine_vis_);
        //output("smoothpoint").set(smoothpoint_vis_);
        output("topology_vis").set(topology_vis_);


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
        for (int j = 0; j < size; j++) {
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

    void LoadReferenceBreaklineNode::process() {
        auto filepath = param<std::string>("filepath");

        PointCollection Vertices_;
        LineString polyline_;

        std::cout << "start load Reference Breakline" << std::endl;
        std::ifstream infile;
        infile.open(filepath);
        if (!infile) {
            std::cout << "LoadReferenceBreaklineNode:: filepath error" << std::endl;
            return;
        }

        std::string xyz;
        while (!infile.eof()) {
            std::getline(infile, xyz);
            auto elems = masb::split2nums(xyz, ',');
            float x = ::atof(elems[0].c_str());
            float y = ::atof(elems[1].c_str());
            float z = ::atof(elems[2].c_str());
            Vertices_.push_back({ x,y,z });
            polyline_.push_back({ x,y,z });
        }
        infile.close();

        std::cout << "Reference breakline has " << Vertices_.size() << " points" << std::endl;

        output("ReferenceVertices").set(Vertices_);
        output("ReferenceBreakline").set(polyline_);
    }

    void LoadTruePositiveVerticesNode::process() {
        auto filepath = param<std::string>("filepath");

        PointCollection Vertices_;
        LineString polyline_;

        std::cout << "start load TruePositiveVertices" << std::endl;
        std::ifstream infile;
        infile.open(filepath);
        if (!infile) {
            std::cout << "LoadTruePositiveVerticesNode:: filepath error" << std::endl;
            return;
        }

        std::string dummyLine;
        int size;
        int i = 0;
        while (i < 1) {
            std::getline(infile, dummyLine);
            i++;
        }

        std::string xyz;
        while (!infile.eof()) {
            std::getline(infile, xyz);
            auto elems = masb::split2nums(xyz, ',');
            float x = ::atof(elems[0].c_str());
            float y = ::atof(elems[1].c_str());
            float z = ::atof(elems[2].c_str());
            Vertices_.push_back({ x,y,z });
            polyline_.push_back({ x,y,z });
        }
        infile.close();

        std::cout << "Loaded breakline has " << Vertices_.size() << " points" << std::endl;

        output("ExtractedTruePositiveVertices").set(Vertices_);
        output("ExtractedTruePositiveBreakline").set(polyline_);
    }

    void SelectTestBreaklineNode::process() {
        auto BreaklineID = param<int>("BreaklineID");
        auto ExtractedBreakline = input("ExtractedBreakline").get<LineStringCollection>();

        if (BreaklineID >= ExtractedBreakline.size()) {
            std::cout << "conpareID error, the maximum id is" << ExtractedBreakline.size() << std::endl;
            return;
        }
        auto TestBreakline_ = ExtractedBreakline[BreaklineID];
        std::cout << "The test breaklline ID is" << BreaklineID
            << "has " << TestBreakline_.size() << " points\n"
            << "the x;y;z coordinate is:\n";
        for (auto &p : TestBreakline_) {
            std::cout << p[0] << ";" << p[1] << ";" << p[2] << "\n";
        }
        output("TestBreakline").set(TestBreakline_);
    }
    void BreaklineValidationNode::process() {
        auto ReferenceVertices = input("ReferenceVertices").get<PointCollection>();
        auto TPVertices_ = input("TruePositiveVertices").get<PointCollection>();

        masb::PointList ReferenceLine;
        ReferenceLine.reserve(ReferenceVertices.size());
        for (auto&v : ReferenceVertices)
            ReferenceLine.push_back(masb::Point(v.data()));

        masb::PointList TestLine;
        for (auto &p : TPVertices_) {
            TestLine.push_back(masb::Point(p.data()));
        }
        ridge::breaklineValidateProcess(TestLine, ReferenceLine);
    }
}
