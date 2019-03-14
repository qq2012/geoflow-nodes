#include "Ma_geometry_processing.hpp"

#include <kdtree2/kdtree2.hpp>

#include <cmath>
#include <iostream>
namespace masb {
    maGeom_result Geom4pt(Vector &p_norm, Vector &q_norm) {
        maGeom_result re;
        //build in cross product function ????????????????????
        /*
        Scalar x = p_norm[1] * q_norm[2] - q_norm[1] * p_norm[2];
        Scalar y = -(p_norm[0] * q_norm[2] - q_norm[0] * p_norm[2]);
        Scalar z = p_norm[0] * q_norm[1] - q_norm[0] * p_norm[1];
        re.bisector = { x,y,z };
        */
        re.bisector = Vrui::Geometry::cross(p_norm, q_norm);
        re.bisector.normalize();//??????2 orientation???
        float cos_SepAng = p_norm * q_norm / p_norm.abs() / q_norm.abs();
        re.SepAng = acos(cos_SepAng);//????value range???
        return re;
    }

    void compute_ma_geometry(ma_data &madata, ma_Geometry &maGeometry) {
        Vector p_norm, q_norm;

        #pragma omp parallel for private(p_norm,q_norm)
        for (int i = 0; i < madata.m*2; ++i) {
            if (madata.ma_qidx[i] != -1) {
                if (i<madata.m)
                    p_norm = (*madata.normals)[i];
                else
                    p_norm = -(*madata.normals)[i - madata.m];
                q_norm = (*madata.normals)[madata.ma_qidx[i]];
                maGeom_result res = Geom4pt(p_norm, q_norm);
                maGeometry.ma_bisector[i] = res.bisector;
                maGeometry.ma_SeperationAng[i] = res.SepAng;
            }
            else {
                //std::cout << "There is a all not shriking\n";
            }
        } 
    }

}
