
#include "Ma_Segmentation_processing.h"
#include <math.h>

using namespace masb;
using namespace std;


//static const double PI = 3.14159265358979323846264338327950288;

void MaSeg_power::update() {
    this->bisec_thres = cos(this->bisecavg_thres / 180.0)*PI;
    this->bisecavg_thres = cos(this->bisecavg_thres / 180.0)*PI;
    this->bisecdiff_thres = cos(this->bisecdiff_thres / 180.0)*PI;
    this->theta_thres = this->theta_thres / 180.0 *PI;
    this->spokecross_thres = cos(this->spokecross_thres / 180.0)*PI;
}

void MaSegProcess::processing(ma_data &madata,ma_Geometry &maGeometry) {
    MaSegProcess SegProcesser;
    SegProcesser.power.update();

}
