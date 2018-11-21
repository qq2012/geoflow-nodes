#include "imgui.h"
#include "geoflow.hpp"

#include <compute_ma_processing.h>
#include <compute_normals_processing.h>

using namespace geoflow;

class ComputeMedialAxisNode:public Node {
  public:
  masb::ma_parameters params;
  float interval = 2;
  double zero=0,pi=3.14;

  ComputeMedialAxisNode(NodeManager& manager):Node(manager) {
    add_input("points", TT_point_collection);
    add_input("normals", TT_vec3f);
    add_output("ma_coords", TT_point_collection);
    add_output("ma_qidx", TT_vec1i);
    add_output("ma_is_interior", TT_vec1i);
  }
  void gui(){
    ImGui::SliderFloat("initial_radius", &params.initial_radius, 0, 1000);
    ImGui::SliderScalar("denoise_preserve", ImGuiDataType_Double, &params.denoise_preserve, &zero, &pi);
    ImGui::SliderScalar("denoise_planar", ImGuiDataType_Double, &params.denoise_planar, &zero, &pi);
    ImGui::Checkbox("nan_for_initr", &params.nan_for_initr);
  }
  void process();
};

class ComputeNormalsNode:public Node {
  public:
  masb::normals_parameters params;
  float interval = 2;

  ComputeNormalsNode(NodeManager& manager):Node(manager) {
    add_input("points", TT_point_collection);
    add_output("normals", TT_vec3f);
  }
  void gui(){
    ImGui::SliderInt("K", &params.k, 1, 100);
  }
  void process();
};