#include "polylineSmooth.h"

#include <cmath>
#include <iostream>
using namespace ridge;

inline void printPolyline(const masb::PointList &polyline) {
    std::cout << "polyline coords: ";
    for (auto &p : polyline) {
        std::cout << "(" << p[0] << "," << p[1] << "," << p[2] << "), ";
    }
    std::cout << "\n";// std::endl;
}
inline float AngleCalaulater(const masb::Point &p1, const masb::Point &c, const masb::Point &p2) {
    auto cp1 = p1 - c;
    auto cp2 = p2 - c;
    cp1.normalize();
    cp2.normalize();
    float cosAng = cp1 * cp2;
    if (cosAng > 1.0 || cosAng < -1.0) {
        std::cout << "cosAng error"
            <<"cosAng"<< cosAng
            << " cp1(" << cp1[0] << "," << cp1[1] << "," << cp1[2] << ")"
            << " cp2(" << cp2[0] << "," << cp2[1] << "," << cp2[2] << ")"
            << std::endl;
    }
    float Ang = acos(cosAng);//value range: (0 - PI rad)
    return Ang;
}

PolyineList ridge::polylineSmooth(const PolyineList &polylines, const float sharpAnlge) {
    //std::cout << "test AngleCalaulater()" << std::endl;
    //AngleCalaulater(masb::Point({ 0,0,1 }), masb::Point({ 0,0,0 }), masb::Point({ 1,0,0 }));

    PolyineList smoothedPolylines = polylines;

    for (auto &a_polyline : smoothedPolylines) {
        std::vector<float> AngList;
        
        std::cout << "a_polyline size" << a_polyline.size() << std::endl;

        for (int idx = 1; idx < a_polyline.size()-1; ++idx) {
            masb::Point p1 = a_polyline[idx - 1];
            masb::Point c = a_polyline[idx];
            masb::Point p2 = a_polyline[idx + 1];
            auto Ang = AngleCalaulater(p1, c, p2);
            AngList.push_back(Ang);
        }
        int minElementIndex = std::min_element(AngList.begin(), AngList.end()) - AngList.begin();
        float minAng = AngList[minElementIndex];
        //std::cout << "minAng degree " << minAng / PI * 180 
        //    <<" minElementIndex"<< minElementIndex << std::endl;
        while (minAng < sharpAnlge){            
            int pt_Idx = minElementIndex + 1;
            int prePt_Idx = pt_Idx - 1;
            int aftPt_Idx = pt_Idx + 1;
            masb::Point prePt = a_polyline[prePt_Idx];
            masb::Point pt = a_polyline[pt_Idx];
            masb::Point aftPt = a_polyline[aftPt_Idx];
            masb::Point newPt1((prePt[0] + pt[0]) / 2.0, (prePt[1] + pt[1]) / 2.0, (prePt[2] + pt[2]) / 2.0);
            masb::Point newPt2((pt[0] + aftPt[0]) / 2.0, (pt[1] + aftPt[1]) / 2.0, (pt[2] + aftPt[2]) / 2.0);
            //update polyline
            a_polyline[pt_Idx] = newPt1;
            a_polyline.insert(a_polyline.begin() + aftPt_Idx, newPt2);
            //update angleList
            auto ang1 = AngleCalaulater(prePt, newPt1, newPt2);
            auto ang2 = AngleCalaulater(newPt1, newPt2, aftPt);
            AngList[minElementIndex] = ang1;
            AngList.insert(AngList.begin() + minElementIndex + 1, ang2);
            //update min Ang
            minElementIndex = std::min_element(AngList.begin(), AngList.end()) - AngList.begin();
            minAng = AngList[minElementIndex];
            //printPolyline(a_polyline);
            //std::cout << "minAng degree " << minAng / PI * 180 
            //    << " minElementIndex" << minElementIndex << std::endl;
        }
        std::cout << "a_polyline size" << a_polyline.size() << std::endl;
    }
    return smoothedPolylines;
}