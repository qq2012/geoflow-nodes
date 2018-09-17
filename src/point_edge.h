#pragma once

#include <iostream>
#include <fstream>

#include <LASlib/lasreader.hpp>
#include <LASlib/laswriter.hpp>

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/linestring.hpp>
#include <boost/geometry/geometries/polygon.hpp>

// CGAL
#include <CGAL/number_utils.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/pca_estimate_normals.h>
#include <CGAL/mst_orient_normals.h>
#include <CGAL/property_map.h>
#include <CGAL/Shape_detection_3.h>
#include <CGAL/regularize_planes.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Search_traits_3.h>
#include <CGAL/linear_least_squares_fitting_3.h>
#include <CGAL/Arr_linear_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Arr_extended_dcel.h>
#include <CGAL/Arr_observer.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Exact_rational.h>
#include <CGAL/Cartesian_converter.h>
#include <utility> // defines std::pair

// line simplification
#include <CGAL/Polyline_simplification_2/simplify.h>

// #include "line_shape.cpp"
#include "region_growing.h"
// Types
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_3 Point;
typedef Kernel::Vector_3 Vector;
typedef Kernel::Vector_2 Vector_2;
typedef Kernel::Plane_3 Plane;
typedef Kernel::Line_3 Line;
// Point with normal vector stored in a std::pair.
typedef CGAL::cpp11::tuple<Point, Vector, int, bool, double, int, bool, double, int> PNL;
typedef CGAL::Nth_of_tuple_property_map<0, PNL> Point_map;
typedef CGAL::Nth_of_tuple_property_map<1, PNL> Normal_map;
typedef CGAL::Nth_of_tuple_property_map<2, PNL> Label_map;
typedef CGAL::Nth_of_tuple_property_map<3, PNL> IsWall_map;
typedef CGAL::Nth_of_tuple_property_map<4, PNL> LineFit_map;
typedef CGAL::Nth_of_tuple_property_map<5, PNL> JumpCount_map;
typedef CGAL::Nth_of_tuple_property_map<6, PNL> IsStep_map;
typedef CGAL::Nth_of_tuple_property_map<7, PNL> JumpEle_map;
typedef CGAL::Nth_of_tuple_property_map<8, PNL> Id_map;
typedef std::vector<PNL>                        PNL_vector;
// Concurrency
#ifdef CGAL_LINKED_WITH_TBB
typedef CGAL::Parallel_tag Concurrency_tag;
#else
typedef CGAL::Sequential_tag Concurrency_tag;
#endif
// In Shape_detection_traits the basic types, i.e., Point and Vector types
// as well as iterator type and property maps, are defined.
typedef CGAL::Shape_detection_3::Shape_detection_traits
  <Kernel, PNL_vector, Point_map, Normal_map>      Traits;
typedef CGAL::Shape_detection_3::Efficient_RANSAC<Traits>    Efficient_ransac;
typedef CGAL::Shape_detection_3::Region_growing<Traits>      Region_growing;
typedef CGAL::Shape_detection_3::Plane<Traits>               SCPlane;

// search tree
// typedef boost::tuple<Point_3,int>                           Point_and_int;
// typedef CGAL::Random_points_in_cube_3<Point_3>              Random_points_iterator;
typedef CGAL::Search_traits_3<Kernel>                       Traits_base;
typedef CGAL::Search_traits_adapter<PNL,
  CGAL::Nth_of_tuple_property_map<0, PNL>,
  Traits_base>                                              TreeTraits;
// typedef CGAL::Search_traits_3<SCK> TreeTraits;
typedef CGAL::Orthogonal_k_neighbor_search<TreeTraits> Neighbor_search;
typedef Neighbor_search::Tree Tree;

// least squares stuff
typedef CGAL::Simple_cartesian<double> SCK;
typedef SCK::Point_3 Point_SCK;
typedef SCK::Line_3 Line_SCK;


// 2D arrangement
// typedef CGAL::Exact_rational Number_type;
// typedef double Number_type;
// typedef CGAL::Cartesian<Number_type>           AK;
typedef CGAL::Exact_predicates_exact_constructions_kernel   AK;
typedef CGAL::Arr_linear_traits_2<AK>                 Traits_2;
typedef Traits_2::Point_2                             Point_2;
typedef CGAL::Polygon_2<AK>                           Polygon_2;
// typedef Kernel::Point_2 PolygonPoint;
typedef Traits_2::Segment_2                           Segment_2;
typedef Traits_2::Ray_2                               Ray_2;
typedef Traits_2::Line_2                              Line_2;
typedef Traits_2::X_monotone_curve_2                  X_monotone_curve_2;
struct FaceInfo {
  bool is_finite=false;
  float elevation_avg=0;
};
typedef CGAL::Arr_face_extended_dcel<Traits_2, FaceInfo>   Dcel;
typedef CGAL::Arrangement_2<Traits_2, Dcel>           Arrangement_2;
typedef Arrangement_2::Vertex_handle                  Vertex_handle;
typedef Arrangement_2::Halfedge_handle                Halfedge_handle;
typedef Arrangement_2::Face_handle                    Face_handle;
typedef Arrangement_2::Vertex_const_handle            Vertex_const_handle;
typedef Arrangement_2::Halfedge_const_handle          Halfedge_const_handle;
typedef Arrangement_2::Face_const_handle              Face_const_handle;
typedef CGAL::Arr_accessor<Arrangement_2>             Arr_accessor;
typedef Arr_accessor::Dcel_vertex                     DVertex;
typedef Arrangement_2::Face                           Face;
// typedef CGAL::Cartesian_converter<IK,EK>                         IK_to_EK;
// typedef CGAL::Cartesian_converter<Number_type,SCK>          ToInexact;

// inline Number_type arr_s(double v){return static_cast<Number_type>(100*v);}

// An arrangement observer, used to receive notifications of face splits and
// to update the indices of the newly created faces.
class Face_index_observer : public CGAL::Arr_observer<Arrangement_2>
{
private:
  int     n_faces;          // The current number of faces.
public:
  Face_index_observer (Arrangement_2& arr) :
    CGAL::Arr_observer<Arrangement_2> (arr),
    n_faces (0)
  {
    CGAL_precondition (arr.is_empty());
    arr.unbounded_face()->data().is_finite=false;
    n_faces++;
  }
  virtual void after_split_face (Face_handle old_face,
                                 Face_handle new_face, bool )
  {
    // Assign index to the new face.
    if(n_faces == 1)
      new_face->data().is_finite = true;
    else if(old_face->data().is_finite)
      new_face->data().is_finite = true;
    else
      new_face->data().is_finite = false;
    n_faces++;
  }
};

// this is for the cgal regularize_plane function
class Index_map
{
  boost::shared_ptr<std::vector<int> > m_indices;
  // const std::vector<int> m_indices;
public:
  typedef int value_type; ///< Index of the shape (-1 if the point is not assigned to any shape).
  typedef int reference;
  typedef boost::readable_property_map_tag category;

  Index_map(){}

  /*!
    Constructs a property map to map points to their associated shape.
    \note `shapes` must be a range of shapes detected using `points`.
    \tparam ShapeRange an `Iterator_range` with a bidirectional
    constant iterator type with value type
    `boost::shared_ptr<CGAL::Shape_detection_3::Shape_base<Traits> >`.
    */
  Index_map (const std::vector<int> intvector)
    : m_indices (new std::vector<int>(intvector))
  {
  }

  inline friend int get (const Index_map& pm, const size_t& k)
  {
    return (*(pm.m_indices))[k];
  }

};

namespace bg = boost::geometry;
typedef bg::model::d2::point_xy<double> point_type;
typedef bg::model::point<double, 3, bg::cs::cartesian> point_type_3d;
typedef bg::model::segment<point_type> segment;

struct config {
  public:
  int metrics_plane_min_points = 40;
  float metrics_plane_epsilon = 0.2;
  float metrics_plane_normal_threshold = 0.95;
  float metrics_is_wall_threshold = 0.1;
  int metrics_k_linefit = 30;
  int metrics_k_jumpcnt_elediff = 10;

  int classify_jump_count_min = 1;
  int classify_jump_count_max = 5;
  float classify_line_dist = 0.005;
  float classify_jump_ele = 1.0;

  float linedetect_dist_threshold = 0.3;
  int linedetect_min_segment_count = 20;
  int linedetect_k = 20;
};

typedef std::vector<std::array<float,2>> vec2f;

void pc_in_footprint(std::string las_filename, std::vector<bg::model::polygon<point_type>> &footprint, std::vector<PNL_vector> &points_vec) ;
void compute_metrics(PNL_vector &points, config = config()) ;
void classify_edgepoints(std::vector<linedect::Point> &edge_points, PNL_vector &points, config = config()) ;
void detect_lines(std::vector<std::pair<Point,Point>> & edge_segments, std::vector<linedect::Point> &edge_points, config = config()) ;
void build_arrangement(bg::model::polygon<point_type> &footprint, std::vector<std::pair<Point,Point>> & edge_segments, Arrangement_2 &arr);
void process_arrangement(PNL_vector& points, Arrangement_2& arr);
void arrangementface_to_polygon(Face_handle face, vec2f& polygons);