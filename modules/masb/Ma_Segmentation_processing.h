
#ifndef MA_SEGMENTATION_PROCESSING_
#define MA_SEGMENTATION_PROCESSING_

#include "madata.h"

#include <CGAL/property_map.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Search_traits_3.h>
#include <CGAL/linear_least_squares_fitting_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Search_traits_adapter.h>
#include <CGAL/Plane_3.h>



namespace masb{
    
  typedef CGAL::Exact_predicates_inexact_constructions_kernel cgal_kernel;
  typedef cgal_kernel::Point_3 cgalPoint;
  //typedef cgal_kernel::Vector_3 Vector;
  //typedef cgal_kernel::Plane_3 Plane;
  
  using namespace std;
  class MaSeg_power{
  public:
    //maSeg_power
    float bisec_thres = 10.0;
    float bisecavg_thres = 2.0;
    float bisecdiff_thres = 5.0;
    float theta_thres = 10.0;
    float secspokecnt_thres = 10;
    float balloverlap_thres = 10;
    int k = 10;
    bool only_interior = true;
    std::string method = "bisec";
    int mincount = 10;
    int maxcount = 1000;
    float spokecross_thres = 5.0;
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
      
    typedef std::pair<cgalPoint,size_t> point_index;
    /*
    typedef CGAL::Search_traits_3<cgal_kernel>  Traits_base;
    typedef CGAL::Search_traits_adapter<point_index,
                                      CGAL::First_of_pair_property_map<point_index>,
                                      Traits_base>  TreeTraits;
    typedef CGAL::Orthogonal_k_neighbor_search<TreeTraits> Neighbor_search;
    typedef Neighbor_search::Tree Tree;
    */
    /*
    vector<Vector> normals;
    Tree tree;
    vector<bool> point_seed_flags;
    size_t region_counter=1;
    */
  public:

    vector<point_index> indexed_points;
    //MaSeg_result
    vector<size_t> point_segment_idx; // 0=unsegmented, maybe put this on the heap...
    //unordered_map<size_t, Plane> segment_shapes;
    MaSeg_power power;
  
    //PlaneDetector(vector<Point> &points, vector<Vector> &normals);
    //vector<size_t> get_point_indices(size_t shape_id);
    void processing(ma_data &madata, ma_Geometry &maGeometry);

  private:
    //inline Plane fit_plane(Neighbor_search search_result);
    inline bool valid_candidate_bisec();
    vector<size_t> findNei(size_t seed_idx);
    void grow_region(size_t seed_idx);

 };
}

#endif
