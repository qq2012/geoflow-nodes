#ifndef EXTRACTCANDIDATEPT_
#define EXTRACTCANDIDATEPT_

#include <iostream>
#include "madata.h"
#include "Ma_utility.h"

//using namespace masb;
namespace masb {
    /*
    struct ExtractCandidatePt_pram {
        std::string method;

    };
    
    class ExtractCandidatePt {
    public:
        int i;
        void processing();
    };
    */
    class ExtractCandidatePt_pram {
    public:
        std::string method = "fixe radius";
        int k_neighbours;
        float SearchRadius = 45;
        float deviationAng_thres = 15.0;
        ExtractCandidatePt_pram() {
            this->deviationAng_thres = cos((this->deviationAng_thres / 180.0)*PI);
        }
    };
    class ExtractCandidatePt {
    public:
        ExtractCandidatePt_pram power;

        PointList can_pt_r;
        PointList can_pt_cos;
        PointList can_pt_bisector_avg;
        intList seg_id;
        VectorList direction;
        VectorList cp;
        VectorList cq;
        //std::vector<PointList> candidate_r, candidate_cos;
        //std::vector <VectorList> candidate_dir;
        //size_t candidate_size = 0;

        /*
            auto pointCloud_ptcollection = input("pointCloud").get<PointCollection>();
    auto candidate_r_ptcollection = input("candidate").get<PointCollection>();
    auto directon_vec3f = input("directon").get<vec3f>();
    auto seg_id_vec1i = input("seg_id").get<vec1i>();
    auto bisector_p_vec3f = input("bisector_p").get<vec3f>();
    auto bisector_q_vec3f = input("bisector_q").get<vec3f>();
    auto adjacency = input("adjacency").get<ridge::int_pair_vec>();
        */


        void processing(ExtractCandidatePt_pram & power, mat_data &madata, ma_Geometry &maGeometry, Sheet_idx_List &sheets,
            PointCloud &PointCloud);
    private:
        size_tList find_trace(ExtractCandidatePt_pram &power, size_t &pt_idx,
            mat_data &maData, ma_Geometry &maGeometry, intList &visitFlag);
        inline bool validateCandidate(float deviationAng_thres, Vector &vec1, Vector &vec2);
    };
}
#endif