#include "masb_nodes.hpp"

namespace geoflow::nodes::mat {   
  NodeRegisterHandle create_register() {
    auto R = NodeRegister::create("MAT");
    R->register_node<ComputeMedialAxisNode>("ComputeMedialAxisNode");
    R->register_node<ComputeNormalsNode>("ComputeNormalsNode");

    R->register_node<PLYLoaderNode>("PLYLoaderNode");
    R->register_node<MaGeometryNode>("MaGeometryNode");
    R->register_node<FilterRNode>("FilterRNode");
    R->register_node<MedialSegmentNode>("MedialSegmentNode");
    R->register_node<MaPt_in_oneTraceNode>("MaPt_in_oneTraceNode");
    R->register_node<ExtractCandidatePtNode>("ExtractCandidatePtNode");
    //R->register_node<ReadCandidatePtNode>("ReadCandidatePtNode");
    //R->register_node<ReadCandidatePtWithBisecNode>("ReadCandidatePtWithBisecNode");
    R->register_node<ConnectCandidatePtNode>("ConnectCandidatePtNode");
    R->register_node<ConnectCandidatePtPolyfitNode>("ConnectCandidatePtPolyfitNode");
    //R->register_node<ReadSegmentRestltNode>("ReadSegmentRestltNode");
    //R->register_node<ReadJunctionPtNode>("ReadJunctionPtNode");
    //R->register_node<ReadAdjacencyNode>("ReadAdjacencyNode");
    R->register_node<LoadReferenceBreaklineNode>("LoadReferenceBreaklineNode");
    R->register_node<BreaklineValidationNode>("BreaklineValidationNode");
    R->register_node<SelectTestBreaklineNode>("SelectTestBreaklineNode");
    R->register_node<LoadTruePositiveVerticesNode>("LoadTruePositiveVerticesNode");

    R->register_node<adjacencyNode>("adjacencyNode");
    R->register_node<PolylineSmothNode>("PolylineSmothNode");
    R->register_node<ExtractCandidatePtAllAtomsNode>("ExtractCandidatePtAllAtomsNode");
    R->register_node<PolyLines3D2objNode>("PolyLines3D2objNode");
    R->register_node<PolylineBSplineSmothNode>("PolylineBSplineSmothNode");
    R->register_node<vectorVisNode>("vectorVisNode");
    //R->register_node<SegmentMakerNode>("SegmentMaker");
    //R->register_node<TestPointsNode>("TestPoints");
    //R->register_node<RegionGrowMedialAxisNode>("RegionGrowMedialAxis");

    return R;
  }

}