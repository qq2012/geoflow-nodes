
#ifndef MA_SEGMENTATION_PROCESSING_
#define MA_SEGMENTATION_PROCESSING_

#include "madata.h"


namespace masb{
  enum METHOD { bisector, radius, thirdopt };
  using namespace std;
  class MaSeg_power{
  public:
    float bisec_thres = 10.0;
    float bisecavg_thres = 2.0;
    float bisecdiff_thres = 5.0;
    float theta_thres = 10.0;
    float secspokecnt_thres = 10;
    float balloverlap_thres = 10;
    int k_neib = 15;
    bool only_interior = true;
    METHOD method = bisector;
    int mincount = 10;
    int maxcount = 1000;
    float spokecross_thres = 5.0;
    float seed_radius_thres = 2.0;
    // mask = None
    void update();
  };

  class MaSegProcess{
      //typedef vector<size_t> sheet;
  public:
      //MaSeg_result
      vector<long long int> point_segment_idx; // 0=unsegmented, maybe put this on the heap...
      Sheet_idx_List shape;//does not contain 0
      vector<int> shape_inout;
      void processing(MaSeg_power &power, mat_data &madata, ma_Geometry &maGeometry);

  private:
      size_t sheet_counter = 1;
      size_t size;
      inline size_t findseed8r(float seed_radius_thres,floatList *ma_radius);
      inline size_t findseed(size_t initial_seed_idx);
      inline bool if_all_segmented();
      inline bool valid_candidate_bisec(float bisec_thres,size_t idx1, size_t idx2, ma_Geometry &maGeometry);
      bool validateCandidate(MaSeg_power &power,size_t idx1, size_t idx2, mat_data &madata, ma_Geometry &maGeometry);
      void grow_sheet(MaSeg_power &power,size_t initial_seed_idx, mat_data &madata, ma_Geometry &maGeometry);
 };
}

#endif
