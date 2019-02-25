
#ifndef MA_SEGMENTATION_PROCESSING_
#define MA_SEGMENTATION_PROCESSING_

#include "madata.h"


namespace masb{
     
  using namespace std;
  class MaSeg_power{
  public:
    float bisec_thres = 10.0;
    float bisecavg_thres = 2.0;
    float bisecdiff_thres = 5.0;
    float theta_thres = 10.0;
    float secspokecnt_thres = 10;
    float balloverlap_thres = 10;
    int k_neib = 10;
    bool only_interior = true;
    std::string method = "bisec";
    int mincount = 10;
    int maxcount = 1000;
    float spokecross_thres = 5.0;
    float seed_radius_thres = 2.0;
    // mask = None

    void update();
    /*
      self.p_bisecthres = math.cos((self.p['bisec_thres'] / 180.0) * math.pi)
      self.p_bisecavgthres = math.cos((self.p['bisecavg_thres'] / 180.0) * math.pi)
      self.p_bisecdiffthres = math.cos((self.p['bisecdiff_thres'] / 180.0) * math.pi)
      self.p_normalthres = math.cos((5.0 / 180.0) * math.pi)
      self.p_thetathres_1 = (self.p['theta_thres'] / 180.0) * math.pi # during bisect growing
      self.p_thetathres_2 = (self.p['theta_thres'] / 180.0) * math.pi # during theta growing
      self.p_k = self.p['k']
      self.p_balloverlap_thres = self.p['balloverlap_thres']
      self.p_mincount = self.p['mincount']
      self.p_spokecross_thres = math.cos((self.p['spokecross_thres'] / 180.0) * math.pi)
   */
  };

  class MaSegProcess{
  public:
      MaSeg_power power;
      //MaSeg_result
      vector<size_t> point_segment_idx; // 0=unsegmented, maybe put this on the heap...
      void processing(ma_data &madata, intList &remainingma_in_out, ma_Geometry &maGeometry);

  private:
      size_t sheet_counter = 1;
      size_t size;
      inline size_t findseed8r(floatList *ma_radius);
      inline size_t findseed();
      inline bool if_all_segmented();
      inline bool valid_candidate_bisec(size_t idx1, size_t idx2, ma_Geometry &maGeometry);
      bool validateCandidate(size_t idx1, size_t idx2, ma_data &madata, ma_Geometry &maGeometry);
      void grow_sheet(size_t initial_seed_idx, ma_data &madata, ma_Geometry &maGeometry);
 };
}

#endif
