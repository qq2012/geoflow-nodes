#include "masb_nodes.hpp"

namespace geoflow::nodes::mat {

  NodeRegister create_register() {
    NodeRegister R("MAT");
    R.register_node<ComputeMedialAxisNode>("ComputeMedialAxisNode");
    R.register_node<ComputeNormalsNode>("ComputeNormalsNode");
    R.register_node<MaGeometryNode>("MaGeometryNode");
    R.register_node<FilterRNode>("FilterRNode");
    R.register_node<MedialSegmentNode>("MedialSegmentNode");
    R.register_node<MaPt_in_oneTraceNode>("MaPt_in_oneTraceNode");
    R.register_node<ExtractCandidatePtNode>("ExtractCandidatePtNode");
    R.register_node<ReadCandidatePtNode>("ReadCandidatePtNode");
    R.register_node<ConnectCandidatePtNode>("ConnectCandidatePtNode");
    return R;
  }

}