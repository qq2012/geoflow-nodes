#include <geoflow/core/geoflow.hpp>

#include <compute_ma_processing.h>
#include <compute_normals_processing.h>
#include "Ma_geometry_processing.hpp"
#include "Ma_utility.h"
//#include "Ma_Segmentation_processing.h"


namespace geoflow::nodes::mat {

  class ComputeMedialAxisNode:public Node {
    public:
    masb::ma_parameters params;
    float interval = 2;
    double zero=0,pi=3.14;
    using Node::Node;
    void init() {
      add_input("points", TT_point_collection);
      add_input("normals", TT_vec3f);
      add_output("ma_coords", TT_point_collection);
      add_output("ma_qidx", TT_vec1i);
      add_output("ma_is_interior", TT_vec1i);
      add_output("ma_radius", TT_vec1f);
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
    using Node::Node;
    void init() {
      add_input("points", TT_point_collection);
      add_output("normals", TT_vec3f);
    }
    void gui(){
      ImGui::SliderInt("K", &params.k, 1, 100);
    }
    void process();
  };

  class MaGeometryNode :public Node {
  public:
      using Node::Node;
      void init() {
          add_input("points", TT_point_collection);
          add_input("normals", TT_vec3f);
          //add_input("ma_coords", TT_point_collection);
          add_input("ma_qidx", TT_vec1i);
          //add_input("ma_radius", TT_vec_float);
          //add_input("ma_is_interior", TT_vec1i);
          add_output("ma_SeparationAng", TT_vec1f);
          add_output("ma_bisector", TT_vec3f);
      }
      void process();
  };
  class FilterRNode :public Node {
  public:
      masb::Filter8R params;
      using Node::Node;
      void init() {
          add_input("ma_radius", TT_vec_float);
          add_output("remaining_idx", TT_vec1i);
      }
      void gui() {
          ImGui::SliderFloat("no radius larger than", &params.radius, 0, 1000);
      }
      void process();

  };
  class MedialSegmentNode :public Node {
  public:
      //masb::MaSeg_power params;
      using Node::Node;
      void init() {
          add_input("points", TT_point_collection);
          add_input("ma_coords", TT_point_collection);
          add_input("remaining_idx", TT_vec1i);
          //add_input("ma_qidx", TT_vec1i);
          add_input("ma_is_interior", TT_vec1i);
          add_output("seg_id", TT_vec1i);
      }
      void gui() {
          //ImGui::SliderFloat("no radius larger than", &params.radius, 0, 1000);
      }
      void process();
  };
}