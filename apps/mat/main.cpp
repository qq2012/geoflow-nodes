#include <iostream>
#include <fstream>

#include "imgui.h"
#include <geoflow/gui/flowchart.hpp>

#include <geoflow/core/geoflow.hpp>
#include <las_register.hpp>
#include <cgal_register.hpp>
#include <masb_register.hpp>
#include <array>

namespace gfn = geoflow::nodes;

int main(int ac, const char * av[])
{
    geoflow::NodeManager N;

    // register nodes
    NodeRegister cgal = gfn::cgal::create_register();
    NodeRegister mat = gfn::mat::create_register();
    NodeRegister las = gfn::las::create_register();

    auto las_loader_node = N.create_node(las, "LASLoader", { 100,100 });
    auto normals_node = N.create_node(mat, "ComputeNormalsNode", { 400,0 });
    auto mat_node = N.create_node(mat, "ComputeMedialAxisNode", { 400,200 });
    auto geometry_node = N.create_node(mat, "MaGeometryNode", { 700,50 });
    auto FilterR_node = N.create_node(mat, "FilterRNode", { 700,300 });
    auto Segmen_node = N.create_node(mat, "MedialSegmentNode", { 1000,200 });
    //auto Trace_node = N.create_node(mat, "MaPt_in_oneTraceNode");
    //auto ExtractCandidatePt_node = N.create_node(mat, "ExtractCandidatePtNode");
    //auto ConnectCandidatePt_node = N.create_node(mat, "ConnectCandidatePtNode");


    las_loader_node->set_param("filepath", "C:/Users/wangq/Downloads/thesis/P3/urban_2.las");//"/some/path"
    connect(las_loader_node->output("points"), normals_node->input("points"));
    connect(normals_node->output("normals"), mat_node->input("normals"));
    connect(las_loader_node->output("points"), mat_node->input("points"));
    connect(las_loader_node->output("points"), geometry_node->input("points"));
    connect(normals_node->output("normals"), geometry_node->input("normals"));
    connect(mat_node->output("ma_qidx"), geometry_node->input("ma_qidx"));
    connect(mat_node->output("ma_radius"), FilterR_node->input("ma_radius"));
    connect(geometry_node->output("ma_SeparationAng"), Segmen_node->input("ma_SeparationAng"));
    connect(geometry_node->output("ma_bisector"), Segmen_node->input("ma_bisector"));
    connect(mat_node->output("ma_coords"), Segmen_node->input("ma_coords"));
    connect(mat_node->output("ma_qidx"), Segmen_node->input("ma_qidx"));
    connect(mat_node->output("ma_radius"), Segmen_node->input("ma_radius"));
    connect(FilterR_node->output("remaining_idx"), Segmen_node->input("remaining_idx"));


    launch_flowchart(N, {cgal, las, mat});
}