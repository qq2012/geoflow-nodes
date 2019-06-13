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
        std::cout << "unshrikingGroundPoint size" << madata.m - mat.matSize << std::endl;
        for (auto &p : mat.unshrikingGroundPoint) {
            unshrikingGroundPoint_.push_back({ p[0],p[1],p[2] });
            //std::cout << p[0] << " " << p[1] << " " << p[2]+0.1 << "\n";
        }
            
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
        
        std::cout << "\nMedialSegmentNode::process()" << std::endl;
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
        output("sheets").set(segmentation.shape);
    }
    void adjacencyNode::process(){
        auto mat = input("mat").get<masb::MAT>();
        auto sheets = input("sheets").get<masb::Sheet_idx_List>();

        std::cout << "\nadjacencyNode::process()" << std::endl;
        ridge::adjacency_parameters  params;
        params.searchRadius = param<float>("searchRadius");
        params.adjacency_thresh = param<int>("adjacency_thresh");

        ridge::int_pair_vec adjacency;
        masb::PointList junction;
        ridge::adjacencyProcessing(params, mat, sheets, adjacency, junction);

        vec3f junction_;
        junction_.reserve(junction.size());
        for (auto&p : junction)
            junction_.push_back({ p[0],p[1],p[2] });

        output("sheet-sheet adjacency").set(adjacency);
        output("junction points").set(junction_);
    };
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
        auto mat = input("mat").get<masb::MAT>();
        auto seg_id_vec1i = input("seg_id").get<vec1i>();
        auto pointcloud_PointCollection = input("pointcloud").get<PointCollection>();
        auto unShrinkingPt_PointCollection = input("unShrinking point cloud").get<PointCollection>();
        
        std::cout << "\nExtractCandidatePtNode::process()" << std::endl;

        masb::ExtractCandidatePt_pram  params;
        params.SearchRadius = param<float>("SearchRadius");
        params.deviationAng_thres = cos((param<float>("deviationAng_thres") / 180.0)*PI);
        std::cout << "deviationAng_thres is " << param<float>("deviationAng_thres") << " degree"
            << "cos (deviationAng_thres) is" << params.deviationAng_thres << std::endl;

        params.MaxEdgeBallRadius = param<float>("MaxEdgeBallRadius");
        params.MinEdgeBallRadius = param<float>("MinEdgeBallRadius");
        params.filterDistance = param<float>("filterDistance");
        params.unshrinkingDist = param<float>("unshrinkingDist");
        params.bis_avg_knn = param<int>("bis_avg_knn");

        masb::PointList pointcloud;
        pointcloud.reserve(pointcloud_PointCollection.size());
        for (auto &p : pointcloud_PointCollection)
            pointcloud.push_back(masb::Point(p.data()));

        masb::PointList unShrinkingPt;
        unShrinkingPt.reserve(unShrinkingPt_PointCollection.size());
        for (auto &p : unShrinkingPt_PointCollection)
            unShrinkingPt.push_back(masb::Point(p.data()));

        masb::intList seg_id;
        seg_id.reserve(seg_id_vec1i.size());
        for (auto id : seg_id_vec1i)
            seg_id.push_back(id);

        masb::ExtractCandidatePt extractor;
        extractor.processing(params,mat, pointcloud, seg_id, unShrinkingPt);

        PointCollection edgeBalls_atoms_; 
        vec1i edgeBall_id_;

        PointCollection candidate_r_, candidate_cos_;
        PointCollection candidate_r_bisector_avg_;

        vec1i filter_;
        PointCollection candidate_points_;
        vec1i candidate_points_id_;
        vec1f candidate_radius_;
        vec3f candidate_direction_;

        int s = extractor.edgeBalls.matSize;
        edgeBalls_atoms_.reserve(s);
        for (auto &p : extractor.edgeBalls.atom)
            edgeBalls_atoms_.push_back({ p[0],p[1],p[2] });

        edgeBall_id_.reserve(s);
        for (auto id : extractor.edgeBall_id)
            edgeBall_id_.push_back(id);

        candidate_r_.reserve(s);
        for (auto&p : extractor.can_pt_r)
            candidate_r_.push_back({ p[0],p[1],p[2] });

        candidate_cos_.reserve(s);
        for (auto&p : extractor.can_pt_cos)
            candidate_cos_.push_back({ p[0],p[1],p[2] });

        candidate_r_bisector_avg_.reserve(s);
        for (auto&p : extractor.can_pt_r_bisector_avg)
            candidate_r_bisector_avg_.push_back({ p[0],p[1],p[2] });

        int s2 = extractor.candidate_pt.size();
        filter_.reserve(s2);
        for (auto i : extractor.filter)
            filter_.push_back(i);

        candidate_points_.reserve(s2);
        for (auto &p : extractor.candidate_pt)
            candidate_points_.push_back({ p[0],p[1],p[2] });

        candidate_points_id_.reserve(s2);
        for (auto i : extractor.candidate_id)
            candidate_points_id_.push_back(i);

        candidate_radius_.reserve(s2);
        for (auto r : extractor.candidate_radius)
            candidate_radius_.push_back(r);

        candidate_direction_.reserve(s2);
        for (auto&v : extractor.candidate_direction)
            candidate_direction_.push_back({ v[0],v[1],v[2] });

        std::cout << "ExtractCandidatePtNode:: there are " << s2 << " candidate points" << std::endl;

        output("edgeBallAtom").set(edgeBalls_atoms_);
        output("edgeBall_id").set(edgeBall_id_);
        output("candidate_r").set(candidate_r_);
        output("candidate_cos").set(candidate_cos_);
        output("candidate_r_bisector_avg").set(candidate_r_bisector_avg_);
        output("filter").set(filter_);
        output("candidate_points").set(candidate_points_);
        output("candidate_points_id").set(candidate_points_id_);
        output("candidate_points_radius").set(candidate_radius_);
        output("candidate_points_direction").set(candidate_direction_);
    }
    void ExtractCandidatePtAllAtomsNode::process() {
        auto mat = input("mat").get<masb::MAT>();
        auto seg_id_vec1i = input("seg_id").get<vec1i>();
        auto pointcloud_PointCollection = input("pointcloud").get<PointCollection>();
        auto unShrinkingPt_PointCollection = input("unShrinking point cloud").get<PointCollection>();

        std::cout << "\nExtractCandidatePtNode::process()" << std::endl;

        masb::ExtractCandidatePt_pram  params;
        params.SearchRadius = 0;
        params.deviationAng_thres = 0;
        params.MaxEdgeBallRadius = 0;
        params.MinEdgeBallRadius = 0;
        params.filterDistance = 0;
        params.unshrinkingDist = param<float>("unshrinkingDist");
        params.bis_avg_knn = param<int>("bis_avg_knn");

        masb::PointList pointcloud;
        pointcloud.reserve(pointcloud_PointCollection.size());
        for (auto &p : pointcloud_PointCollection)
            pointcloud.push_back(masb::Point(p.data()));

        masb::PointList unShrinkingPt;
        unShrinkingPt.reserve(unShrinkingPt_PointCollection.size());
        for (auto &p : unShrinkingPt_PointCollection)
            unShrinkingPt.push_back(masb::Point(p.data()));

        masb::intList seg_id;
        seg_id.reserve(seg_id_vec1i.size());
        for (auto id : seg_id_vec1i)
            seg_id.push_back(id);

        masb::PointList candidatePt;
        masb::intList candidatePt_id;
        masb::allAtoms2Candidates(params, mat, pointcloud, seg_id, unShrinkingPt, 
            candidatePt, candidatePt_id);

        PointCollection candidatePt_;
        vec1i candidatePt_id_;

        int s = candidatePt.size();
        candidatePt_.reserve(s);
        for (auto &p : candidatePt)
            candidatePt_.push_back({ p[0],p[1],p[2] });

        candidatePt_id_.reserve(s);
        for (auto id : candidatePt_id)
            candidatePt_id_.push_back(id);

        std::cout << "ExtractCandidatePtAllAtomsNode:: there are " << s << " candidate points" << std::endl;

        output("candidate_points").set(candidatePt_);
        output("candidate_points_id").set(candidatePt_id_);
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
        auto candidate_pt_ptcollection = input("candidate_points").get<PointCollection>();
        auto point_directon_vec3f = input("point_directon").get<vec3f>();
        auto candidate_points_id_vec1i = input("candidate_points_id").get<vec1i>();
        auto adjacency = input("adjacency").get<ridge::int_pair_vec>();
        
        std::cout << "\nConnectCandidatePtNode::process()" << std::endl;

        size_t size = candidate_pt_ptcollection.size();
        masb::PointList candidate_pt;
        candidate_pt.reserve(size);
        for (auto& p : candidate_pt_ptcollection) {
            candidate_pt.push_back(masb::Point(p.data()));
        }
        
        masb::VectorList pt_directon;
        pt_directon.reserve(size);
        for (auto& n : point_directon_vec3f) {
            auto n_normlize = masb::Vector(n.data()).normalize();
            pt_directon.push_back(n_normlize);
        }
        
        masb::intList candidate_pt_id;
        candidate_pt_id.reserve(size);
        for (auto id : candidate_points_id_vec1i) {
            candidate_pt_id.push_back(id);
        }

        ridge::LineSegmentList mstLineSegment;
        masb::intList mstLineSegment_id;
        ridge::PolyineList polylines_maxDistance, polylines_maxAccDist, polylines_maxPtNum;
        masb::intList polyline_id;

        std::cout << "connect candidate points to closest neighbour "
            <<"and symplify as one polyline" << std::endl;
        ridge::connectCandidatePt8MST(candidate_pt, pt_directon, candidate_pt_id,
            mstLineSegment, mstLineSegment_id, polylines_maxDistance, polylines_maxAccDist, 
            polylines_maxPtNum, polyline_id);
 
        ridge::PolyineList linesWithJunction1 = polylines_maxAccDist;
        ridge::FindTopology(linesWithJunction1, polyline_id, adjacency);

        ridge::PolyineList linesWithJunction2 = polylines_maxDistance;
        ridge::FindTopology(linesWithJunction2, polyline_id, adjacency);

        //std::cout << "STARTING METHOD 2 -- NO SEG-ID\n";
        //ridge::connectCandidatePt8MST_nosegid(pointCloud, candidate_r, filter2, segmentList_nosegid);

        LineStringCollection mstLineSegment_;
        mstLineSegment_.reserve(mstLineSegment.size());
        for (auto &pointpair : mstLineSegment) {
            LineString tmp;
            auto p = pointpair.first;
            tmp.push_back({ p[0],p[1],p[2] });
            auto q = pointpair.second;
            tmp.push_back({ q[0],q[1],q[2] });
            mstLineSegment_.push_back(tmp);
        }

        vec1i mstLineSegment_id_;
        mstLineSegment_id_.reserve(mstLineSegment_id.size());
        for (auto i : mstLineSegment_id)
            mstLineSegment_id_.push_back(i);

        LineStringCollection polylines_maxDistance_;
        polylines_maxDistance_.reserve(polylines_maxDistance.size());
        for (auto &a_line : polylines_maxDistance) {
            LineString tmp;
            for (auto &p : a_line) {
                tmp.push_back({ p[0],p[1],p[2] });
            }
            polylines_maxDistance_.push_back(tmp);
        }
        LineStringCollection polylines_maxAccDist_;
        polylines_maxAccDist_.reserve(polylines_maxAccDist.size());
        for (auto &a_line : polylines_maxAccDist) {
            LineString tmp;
            for (auto &p : a_line) {
                tmp.push_back({ p[0],p[1],p[2] });
            }
            polylines_maxAccDist_.push_back(tmp);
        }

        LineStringCollection polylines_maxPtNum_;
        polylines_maxPtNum_.reserve(polylines_maxPtNum.size());
        for (auto &a_line : polylines_maxPtNum) {
            LineString tmp;
            for (auto &p : a_line) {
                tmp.push_back({ p[0],p[1],p[2] });
            }
            polylines_maxPtNum_.push_back(tmp);
        }

        vec1i polyline_id_;
        polyline_id_.reserve(polyline_id.size());
        for (auto i : polyline_id)
            polyline_id_.push_back(i);

        LineStringCollection linesWithJunction_maxAccDist_;
        linesWithJunction_maxAccDist_.reserve(linesWithJunction1.size());
        for (auto &a_line : linesWithJunction1) {
            LineString tmp;
            for (auto &p : a_line) {
                tmp.push_back({ p[0],p[1],p[2] });
            }
            linesWithJunction_maxAccDist_.push_back(tmp);
        }

        LineStringCollection linesWithJunction_maxDistance_;
        linesWithJunction_maxDistance_.reserve(linesWithJunction2.size());
        for (auto&a_line : linesWithJunction2) {
            LineString tmp;
            for (auto &p : a_line) {
                tmp.push_back({ p[0],p[1],p[2] });
            }
            linesWithJunction_maxDistance_.push_back(tmp);
        }


        output("mstLineSegment").set(mstLineSegment_);
        output("mstLineSegment_id").set(mstLineSegment_id_);
        output("polylines_maxDistance").set(polylines_maxDistance_);
        output("polylines_maxAccDist").set(polylines_maxAccDist_);
        output("polylines_maxPtNum").set(polylines_maxPtNum_);
        output("polyline_id").set(polyline_id_);
        output("linesWithJunction_maxAccDist").set(linesWithJunction_maxAccDist_);
        output("linesWithJunction_maxDistance").set(linesWithJunction_maxDistance_);
    }
    void ConnectCandidatePtPolyfitNode::process() {
        auto error_thresh = param<float>("min_error_thresh");

        auto candidate_pt_ptcollection = input("candidate_points").get<PointCollection>();
        auto candidate_points_id_vec1i = input("candidate_points_id").get<vec1i>();
        auto pointcloud_PointCollection = input("pointcloud").get<PointCollection>();
        //auto adjacency = input("adjacency").get<ridge::int_pair_vec>();

        std::cout << "\nConnectCandidatePtPolyfitNode::process()" << std::endl;

        size_t size = candidate_pt_ptcollection.size();
        masb::PointList candidate_pt;
        candidate_pt.reserve(size);
        for (auto& p : candidate_pt_ptcollection) {
            candidate_pt.push_back(masb::Point(p.data()));
        }

        masb::intList candidate_pt_id;
        candidate_pt_id.reserve(size);
        for (auto id : candidate_points_id_vec1i) {
            candidate_pt_id.push_back(id);
        }

        masb::PointList pointcloud;
        pointcloud.reserve(pointcloud_PointCollection.size());
        for (auto &p : pointcloud_PointCollection)
            pointcloud.push_back(masb::Point(p.data()));

        //float error_thresh = 1000000;
        ridge::PolyineList polylines;
        ridge::ConnectCandidatePt8PolynomialFitting(candidate_pt, candidate_pt_id, pointcloud, error_thresh, polylines);

        LineStringCollection polylines_;
        polylines_.reserve(polylines.size());
        for (auto&a_line : polylines) {
            LineString tmp;
            for (auto &p : a_line) {
                tmp.push_back({ p[0],p[1],p[2] });
            }
            polylines_.push_back(tmp);
        }

        output("polylines").set(polylines_);
    }

    void PolylineSmothNode::process() {
        auto polyline_linestringcollection = input("polylines").get<LineStringCollection>();

        float sharpAng_degree = param<float>("sharpAng");
        float sharpAng_rad = sharpAng_degree / 180.0*PI;
        
        ridge::PolyineList polylines;
        polylines.reserve(polyline_linestringcollection.size());
        for (auto &polyline_linestring : polyline_linestringcollection) {
            masb::PointList a_polyline;
            for (auto &pt : polyline_linestring) {
                a_polyline.push_back(masb::Point(pt.data()));
            }
            polylines.push_back(a_polyline);
        }
        /*
        //test
        ridge::PolyineList polylines;
        masb::PointList a_polyline;
        a_polyline.push_back(masb::Point({ 1,2,0 }));
        a_polyline.push_back(masb::Point({ 3,2,0 }));
        a_polyline.push_back(masb::Point({ 5,2,0 }));
        a_polyline.push_back(masb::Point({ 5,0,0 }));
        a_polyline.push_back(masb::Point({ 7,0,0 }));
        a_polyline.push_back(masb::Point({ 9,0,0 }));
        a_polyline.push_back(masb::Point({ 11,0,0 }));
        polylines.push_back(a_polyline);

        LineStringCollection testPolyline_coll;
        LineString testpolyline;
        testpolyline.push_back({ 1,2,0 });
        testpolyline.push_back({ 3,2,0 });
        testpolyline.push_back({ 5,2,0 });
        testpolyline.push_back({ 5,0,0 });
        testpolyline.push_back({ 7,0,0 });
        testpolyline.push_back({ 9,0,0 });
        testpolyline.push_back({ 11,0,0 });
        testPolyline_coll.push_back(testpolyline);
        */
        auto smoothedPolylines = ridge::polylineSmooth(polylines, sharpAng_rad);
        
        LineStringCollection smoothedPolylines_;
        smoothedPolylines_.reserve(smoothedPolylines.size());
        for (auto &a_line : smoothedPolylines) {
            LineString tmp;
            for (auto&p : a_line) {
                tmp.push_back({ p[0],p[1],p[2] });
            }
            smoothedPolylines_.push_back(tmp);
        }
        output("smoothed polyline").set(smoothedPolylines_);
        //output("test polyline").set(testPolyline_coll);
    }

    void PolylineBSplineSmothNode::process() {
        auto polyline_linestringcollection = input("polylines").get<LineStringCollection>();
        
        ridge::PolyineList polylines;
        polylines.reserve(polyline_linestringcollection.size());
        for (auto &polyline_linestring : polyline_linestringcollection) {
            masb::PointList a_polyline;
            for (auto &pt : polyline_linestring) {
                a_polyline.push_back(masb::Point(pt.data()));
            }
            polylines.push_back(a_polyline);
        }
        
        //test
        /*
        ridge::PolyineList polylines;
        masb::PointList a_polyline;
        a_polyline.push_back(masb::Point({ 1,2,13 }));
        a_polyline.push_back(masb::Point({ 3,2,11 }));
        a_polyline.push_back(masb::Point({ 5,2,14 }));
        a_polyline.push_back(masb::Point({ 5,0,25 }));
        a_polyline.push_back(masb::Point({ 7,0,27 }));
        a_polyline.push_back(masb::Point({ 9,0,12 }));
        a_polyline.push_back(masb::Point({ 11,0,15 }));
        polylines.push_back(a_polyline);
        */
        
        //masb::PointList testpolyline;
        //auto smoothedPolylines = ridge::polylineSmooth8EigenSpline(polylines);
        auto smoothedPolylines = ridge::polylineBSplineSmooth(polylines);
        /*
        LineStringCollection testPolyline_radom3d;
        LineString atestPolyline_radom3d;
        for (auto pt : testpolyline)
            atestPolyline_radom3d.push_back({ pt[0],pt[1],pt[2] });
        testPolyline_radom3d.push_back(atestPolyline_radom3d);
        */
        LineStringCollection smoothedPolylines_;
        smoothedPolylines_.reserve(smoothedPolylines.size());
        for (auto &a_line : smoothedPolylines) {
            LineString tmp;
            for (auto&p : a_line) {
                tmp.push_back({ p[0],p[1],p[2] });
            }
            smoothedPolylines_.push_back(tmp);
        }
        output("smoothed polyline").set(smoothedPolylines_);
        //output("test polyline").set(testPolyline_radom3d);
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
    void PolyLines3D2objNode::process() {
        auto polylines_LineStringCollection = input("polylines").get<LineStringCollection>();
        auto filepath = param<std::string>("filepath");
        
        std::ofstream outfile;
        outfile.open(filepath);
        if (!outfile) {
            std::cout << "PolyLines3D2objNode:: filepath error" << std::endl;
            return;
        }

        for (auto &linestring : polylines_LineStringCollection) {
            for (auto&p : linestring) {
                outfile << "v " << p[0]
                    << " " << p[1]
                    << " " << p[2]
                    << "\n";
            }
        }
        int count = 1;
        for (auto& linestring : polylines_LineStringCollection) {
            outfile << "l ";
            int cur_size = linestring.size();
            for (int i = count; i < count + cur_size; ++i) {
                outfile << i << " ";
            }
            outfile << "\n";
            count += cur_size;
        }
        outfile.close();
        //wkt LINESTRING (1 2 0, 4 3 0, 8 9 0)
    }
    void vectorVisNode::process() {
        auto startPoint = input("startPoint").get<PointCollection>();
        auto direction = input("direction").get<vec3f>();

        if (startPoint.size() != direction.size()) {
            std::cout << "vectorVisNode::process() size error" << std::endl;
            return;
        }
        LineStringCollection vis_;
        vis_.reserve(startPoint.size());
        for (int i = 0; i < startPoint.size(); ++i) {
            LineString tmp;
            auto pt = startPoint[i];
            auto v = direction[i];
            tmp.push_back(pt);
            tmp.push_back({ pt[0] + v[0],pt[1] + v[1],pt[2] + v[2] });
            vis_.push_back(tmp);
        }
        output("visulization").set(vis_);
    }
}
