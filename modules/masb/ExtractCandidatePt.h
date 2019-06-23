#ifndef EXTRACTCANDIDATEPT_
#define EXTRACTCANDIDATEPT_

#include <iostream>
#include "madata.h"
#include "Ma_utility.h"

//using namespace masb;
namespace masb {
    struct ExtractCandidatePt_pram{
    public:
        float SearchRadius;// = 45;
        float deviationAng_thres;// cosin value
        float MaxEdgeBallRadius; //filter edge ball with big r
        float MinEdgeBallRadius; //filter edge ball with small r
        float filterDistance; //filter edge ball far away from ground
        float unshrinkingDist;//distance from candidate pt (x,y,z_interpolation) to unshrinking point 
        int bis_avg_knn;
    };

    void allAtoms2Candidates(ExtractCandidatePt_pram & power, MAT &mat, PointList &pointcloud,
        intList &seg_id, PointList &unShrinkingPt, PointList &candidatePt, intList &candidatePt_id);

    class ExtractCandidatePt {
    public:
        MAT edgeBalls;
        intList edgeBall_id;
        PointList can_pt_r;
        PointList can_pt_cos;
        PointList can_pt_r_bisector_avg;

        intList filter;
        PointList candidate_pt;
        intList candidate_id;
        floatList candidate_radius;
        VectorList candidate_direction;

        void processing(ExtractCandidatePt_pram & power, MAT &mat, PointList &pointcloud, 
            intList &seg_id, PointList &unShrinkingPt);
        
    private:
        //PointList remainingPt;
        inline bool validateCandidate(float deviationAng_thres, Vector &vec1, Vector &vec2);
        void EdgeBallDetection(ExtractCandidatePt_pram & power, MAT &mat,int cur_id);
        //void EdgeBall2CandidatePt();
        //void spatialInterpolation(PointList &pointcloud);
        void filterCandidatePt(ExtractCandidatePt_pram & power, PointList &pointcloud, PointList &unShrinkingPt);

    };
}
#endif