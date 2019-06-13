
find_package (Eigen3 3.3 REQUIRED NO_MODULE)
add_library(masb_nodes STATIC
  masb_nodes.cpp
  Ma_geometry_processing.cpp
  Ma_utility.cpp
  Ma_Segmentation_processing.cpp
  ExtractCandidatePt.cpp
  ExtractCandidatePt_AllAtoms.cpp
  MaPt_in_oneTrace.cpp
  ConnectCandidatePt.cpp
  ConnectCandidatePt_PolynomialFitting.cpp
  FindTopology.cpp
  breaklineValidate.cpp
  sheet_adjacency.cpp
  polylineSmooth.cpp
  polylineSmooth_Spline.cpp
  PolynomialFitting.cpp
  Bezier.cpp
  BSpline.cpp
  Curve.cpp
  Vector.cpp
)

include_directories(${EIGEN_INCLUDE_DIR})

target_link_libraries( masb_nodes
  geoflow-core
  masbcpp
  Eigen3::Eigen
)
set_target_properties(
  masb_nodes
  PROPERTIES CXX_STANDARD 17
)
#add_executable (example example.cpp)