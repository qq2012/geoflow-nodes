#include "masb_nodes.hpp"

namespace geoflow::nodes::mat {

  NodeRegister create_register() {
    NodeRegister R("MAT");
    R.register_node<ComputeMedialAxisNode>("ComputeMedialAxisNode");
    R.register_node<ComputeNormalsNode>("ComputeNormalsNode");
    R.register_node<MaGeometryNode>("MaGeometryNode");
    R.register_node<MedialBisecSegmentNode>("MedialBisecSegmentNode");
    return R;
  }

}