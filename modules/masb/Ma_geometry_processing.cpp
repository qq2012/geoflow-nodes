#include "Ma_geometry_processing.hpp"

#include <kdtree2/kdtree2.hpp>

#include <cmath>
#include <iostream>
namespace masb {
    void compute_ma_geometry_new(const ma_data &madata, MAT &mat, const bool isInterior) {
        int startIdx, endIdx;
        if (isInterior) {
            startIdx = 0;
            endIdx = madata.m;
        }
        else {
            startIdx = madata.m;
            endIdx = madata.m * 2;
        }
        for (int i = startIdx; i < endIdx; ++i) {
            int i_ = i % madata.m;
            if (madata.ma_qidx[i] != -1) {
                Point c = (*madata.ma_coords)[i];
                mat.atom.push_back(c);
                mat.radius.push_back((*madata.ma_radius)[i]);
                Point p = (*madata.coords)[i_];
                mat.fp.push_back(p);//feature point p
                Point q = (*madata.coords)[madata.ma_qidx[i]];
                mat.fq.push_back(q);
                if (isInterior)
                    mat.sp_reverse_norm.push_back((*madata.normals)[i_]);
                else
                    mat.sp_reverse_norm.push_back(-(*madata.normals)[i_]);
                auto sp = p - c;
                sp.normalize();
                auto sq = q - c;
                sq.normalize();
                mat.sp.push_back(sp);//spoken vector
                mat.sq.push_back(sq);
                auto b = sp + sq;
                b.normalize();
                mat.bisector.push_back(b);
                float sepAng = std::acos(sp*sq);
                mat.seperationAng.push_back(sepAng);
                auto md = Vrui::Geometry::cross(sp, sq);
                mat.ma_direction.push_back(md);
            }
            else {
                Point u = (*madata.coords)[i_];
                mat.unshrikingGroundPoint.push_back(u);
            }
        }

    }


    maGeom_result Geom4pt(Vector &p_norm, Vector &q_norm) {
        maGeom_result re;
        //build in cross product function ????????????????????
        /*
        Scalar x = p_norm[1] * q_norm[2] - q_norm[1] * p_norm[2];
        Scalar y = -(p_norm[0] * q_norm[2] - q_norm[0] * p_norm[2]);
        Scalar z = p_norm[0] * q_norm[1] - q_norm[0] * p_norm[1];
        re.norm = { x,y,z };
        */
        re.bisector = p_norm + q_norm;
        re.bisector.normalize();
        re.norm = Vrui::Geometry::cross(p_norm, q_norm);
        re.norm.normalize();
        float cos_SepAng = p_norm * q_norm / p_norm.abs() / q_norm.abs();
        re.SepAng = acos(cos_SepAng);//value range: 0-180 degree (0 - PI rad)
        return re;
    }

    void compute_ma_geometry(ma_data &madata, ma_Geometry &maGeometry) {
        //why dose not work!!!!????
        /*
        maGeometry.ma_bisector.reserve(madata.m * 2);
        maGeometry.ma_normal.reserve(madata.m * 2);
        maGeometry.ma_SeperationAng.reserve(madata.m * 2);
        maGeometry.ma_bisector(madata.m * 2);
        maGeometry.ma_normal(madata.m * 2);
        maGeometry.ma_SeperationAng(madata.m * 2);
        */
        Vector cp,cq;
        #pragma omp parallel for private(cp,cq)
        for (int i = 0; i < madata.m * 2; ++i) {
            if (madata.ma_qidx[i] != -1){
                Point p;
                if (i < madata.m){
                    //cp = (*madata.normals)[i];
                    p = (*madata.coords)[i];
                }
                else {
                    //cp = -(*madata.normals)[i - madata.m];
                    p = (*madata.coords)[i - madata.m];
                }
                Point q = (*madata.coords)[madata.ma_qidx[i]];
                Point c = (*madata.ma_coords)[i];
                cq= q - c;
                cq.normalize();
                cp = p - c;
                cp.normalize();
                maGeom_result res = Geom4pt(cp, cq);
                maGeometry.ma_bisector[i] = res.bisector;
                maGeometry.ma_SeperationAng[i] = res.SepAng;
                maGeometry.ma_normal[i] = res.norm;
            }
            else {
                //std::cout << "There is a all not shriking\n";
            }
        }
        maGeometry;
    }

}
