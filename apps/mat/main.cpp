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

    auto las_loader_node = N.create_node(las, "LASLoader", { -300,0 });
    auto mat_node = N.create_node(mat, "ComputeMedialAxisNode", { 0,0 });
    auto normals_node = N.create_node(mat, "ComputeNormalsNode", { -100,0 });

    las_loader_node->set_param("filepath", "/some/path");

    connect(normals_node->output("normals"), mat_node->input("normals"));


    launch_flowchart(N, {cgal, las, mat});
}