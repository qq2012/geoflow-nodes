
#ifndef MA_SEGMENTATION_PROCESSING_
#define MA_SEGMENTATION_PROCESSING_

#include "madata.h"


namespace masb{
  enum METHOD { bisector, spokecross, balloverlap, 
      combinBisecAndSpcros, combinBallAndSpcros};
  using namespace std;
  class MaSeg_power{
  public:
      METHOD method;
      float bisec_thres;// turned into cosin
      float spokecross_thres;// turned into cosin
      float seed_radius_thres;// distance
      float balloverlap_thres;// (r1+r2)/distance12  it is a scalar 
      int mincount;// number
      int maxcount;// number
      int k_neib = 15;//number

      //float bisecavg_thres;// = 2.0;
      //float bisecdiff_thres;// = 5.0;
      //float theta_thres;// = 10.0;
      //float secspokecnt_thres;// = 10;
      //bool only_interior = true;
      // mask = None
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
      vector<long long int> *remaining_idx;
      size_t sheet_counter = 1;
      size_t size;
      inline size_t findseed8r(float seed_radius_thres,floatList *ma_radius);
      inline size_t findseed(size_t initial_seed_idx);
      size_t findseed_random();
      void  remaining_idx_remove_idx(long long int id);
      inline bool if_all_segmented();
      inline bool valid_candidate_bisec(float bisec_thres,size_t idx1, size_t idx2, ma_Geometry &maGeometry);
      inline bool valid_candidate_spokecross(float cosnorm_thres, size_t idx1, size_t idx2, ma_Geometry &maGeometry);
      inline bool valid_candidate_balloverlap(float balloverlap_thres, size_t idx1, size_t  idx2, mat_data &madata);
      bool validateCandidate(MaSeg_power &power,size_t idx1, size_t idx2, mat_data &madata, ma_Geometry &maGeometry);
      void grow_sheet(MaSeg_power &power,size_t initial_seed_idx, mat_data &madata, ma_Geometry &maGeometry);
 };
}

#endif
