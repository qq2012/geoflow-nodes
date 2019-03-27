add_library(masb_nodes STATIC
  masb_nodes.cpp
  Ma_geometry_processing.cpp
  Ma_utility.cpp
  Ma_Segmentation_processing.cpp
  ExtractCandidatePt.cpp
  MaPt_in_oneTrace.cpp
  ConnectCandidatePt.cpp
)
target_link_libraries( masb_nodes
  geoflow-core
  masbcpp
)
set_target_properties(
  masb_nodes
  PROPERTIES CXX_STANDARD 17
)