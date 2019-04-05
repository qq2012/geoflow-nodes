
find_package (Eigen3 3.3 REQUIRED NO_MODULE)
add_library(masb_nodes STATIC
  masb_nodes.cpp
  Ma_geometry_processing.cpp
  Ma_utility.cpp
  Ma_Segmentation_processing.cpp
  ExtractCandidatePt.cpp
  MaPt_in_oneTrace.cpp
  ConnectCandidatePt.cpp
  ConnectCandidatePt_smooth.cpp
  ConnectCandidatePt_spline.cpp
)
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