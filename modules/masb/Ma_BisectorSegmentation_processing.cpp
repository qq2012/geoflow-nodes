/*
#include "Ma_BisectorSegmentation_processing.h"
#include <math.h>

using namespace masb;
using namespace std;


static const double PI = 3.14159265358979323846264338327950288;

void MaBiSeg_power::update() {
    this->bisec_thres = cos(this->bisecavg_thres / 180.0)*PI;
    this->bisecavg_thres = cos(this->bisecavg_thres / 180.0)*PI;
    this->bisecdiff_thres = cos(this->bisecdiff_thres / 180.0)*PI;
    this->theta_thres = this->theta_thres / 180.0 *PI;
    this->spokecross_thres = cos(this->spokecross_thres / 180.0)*PI;
}

void MaBiSegProcess::Seg_process(MaBiSeg_power &input_parameters) {
    MaBiSeg_power power;
    power.update();

}
*/