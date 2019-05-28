
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
  };

  class MaSegProcess{
  public:
      vector<long long int> point_segment_idx; // 0=unsegmented
      Sheet_idx_List shape;//does not contain 0  <<idx1,idx3......>,<idx5,idx6,idx7>......>
      void processing(MaSeg_power &power, MAT &mat);
  private:
      vector<long long int> *remaining_idx;
      size_t sheet_counter = 1;
      size_t size;
      inline size_t findseed(size_t initial_seed_idx);
      //size_t findseed_random();
      //void  remaining_idx_remove_idx(long long int id);
      inline bool if_all_segmented();
      inline bool valid_candidate_bisec(float bisec_thres,size_t idx1, size_t idx2, MAT&mat);
      inline bool valid_candidate_sepAng(float sepAng_thres, size_t idx1, size_t  idx2, MAT&mat);
      inline bool valid_candidate_spokecross(float cosnorm_thres, size_t idx1, size_t idx2, MAT&mat);
      inline bool valid_candidate_balloverlap(float balloverlap_thres, size_t idx1, size_t  idx2, MAT&mat);
      bool validateCandidate(MaSeg_power &power,size_t idx1, size_t idx2, MAT&mat);
      void grow_sheet(MaSeg_power &power,size_t initial_seed_idx, MAT &mat, kdtree2::KDTree *kdtree_ma_atoms);
 };
}

#endif
