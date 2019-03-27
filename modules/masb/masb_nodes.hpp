#include <geoflow/core/geoflow.hpp>

#include <compute_ma_processing.h>
#include <compute_normals_processing.h>
#include "Ma_geometry_processing.hpp"
#include "Ma_utility.h"
#include "Ma_Segmentation_processing.h"
#include "ExtractCandidatePt.h"
#include "MaPt_in_oneTrace.h"
#include "ConnectCandidatePt.h"


namespace geoflow::nodes::mat {

  class ComputeMedialAxisNode:public Node {
    public:
    float interval = 2;
    double zero=0,pi=3.14;
    using Node::Node;
    void init() {
        add_param("initial_radius", (float)200);
        add_param("denoise_preserve", (double)((PI / 180.0) * 20));
        add_param("denoise_planar", (double)((PI / 180.0) * 32));
        add_param("nan_for_initr", (bool)false);
        add_input("points", TT_point_collection);
        add_input("normals", TT_vec3f);
        add_output("ma_coords", TT_point_collection);
        add_output("ma_qidx", TT_vec1i);
        add_output("ma_is_interior", TT_vec1i);
        add_output("ma_radius", TT_vec1f);
        add_output("coords_masb", TT_any);
    }
    void gui(){
      ImGui::SliderFloat("initial_radius", &param<float>("initial_radius"), 0, 1000);
      ImGui::SliderScalar("denoise_preserve", ImGuiDataType_Double, &param<double>("denoise_preserve"), &zero, &pi);
      ImGui::SliderScalar("denoise_planar", ImGuiDataType_Double, &param<double>("denoise_planar"), &zero, &pi);
      ImGui::Checkbox("nan_for_initr", &param<bool>("nan_for_initr"));
    }
    void process();
  };

  class ComputeNormalsNode:public Node {
    public:
    //masb::normals_parameters params;
    float interval = 2;
    using Node::Node;
    void init() {
        add_param("K", (int) 10);
        add_input("points", TT_point_collection);
        add_output("normals", TT_vec3f);
    }
    void gui(){
      //ImGui::SliderInt("K", &params.k, 1, 100);
      ImGui::SliderInt("K", &param<int>("K"), 1, 100);
    }
    void process();
  };

  class MaGeometryNode :public Node {
  public:
      using Node::Node;
      void init() {
          add_input("points", TT_point_collection);
          add_input("normals", TT_vec3f);
          add_input("ma_coords", TT_point_collection);
          add_input("ma_qidx", TT_vec1i);
          add_output("ma_SeparationAng", TT_vec1f);
          add_output("ma_bisector", TT_vec3f);
          add_output("bisec",TT_line_string_collection);
          add_output("ma_normal", TT_vec3f);
      }
      void process();
  };
  class FilterRNode :public Node {
  public:
      //float radius = 200.00;
      using Node::Node;
      void init() {
          add_param("filterRadius", (float) 100.00);
          //add_param("filterRadius", param<float>("initial_radius"));
          add_input("ma_radius", TT_vec1f);
          add_output("remaining_idx", TT_vec1i);
      }
      void gui() {
          //ImGui::SliderFloat("no radius larger than", &radius, 0, 1000);
          //ImGui::SliderFloat("filterRadius", &param<float>("filterRadius"), 0, param<float>("initial_radius"));
          ImGui::SliderFloat("filterRadius", &param<float>("filterRadius"), 0, 1000);
      }
      void process();

  };
  class MedialSegmentNode :public Node {
  private:
      masb::METHOD current_method = masb::bisector;
  public:
      using Node::Node;
      void init() {
          add_param("mincount", (int)10);
          add_param("maxcount", (int)1000);
          //add_param("method", "bisec");//why not std::string?????
          add_param("bisec_thres", (float)10.0);
          add_param("balloverlap_thres", (float)10);

          add_input("remaining_idx", TT_vec1i);
          add_input("ma_coords", TT_point_collection);
          add_input("ma_qidx", TT_vec1i);
          add_input("ma_radius", TT_vec1f);
          add_input("ma_SeparationAng", TT_vec1f);
          add_input("ma_bisector", TT_vec3f);
          add_input("ma_normal", TT_vec3f);

          add_output("seg_id", TT_vec1i);
          add_output("sheet_all", TT_any);//#########  todo type -- line string clooection
          add_output("ma_coords", TT_point_collection);
          add_output("madata_in",TT_any);
          add_output("madata_out", TT_any);
          add_output("maGeometry_in",TT_any);
          add_output("maGeometry_out", TT_any);
          add_output("seg_in", TT_any);
          add_output("seg_out", TT_any);
          add_output("sheet_in", TT_any);
          add_output("sheet_out",TT_any);
      }
      void gui() {
          ImGui::SliderInt("mincount", &param<int>("mincount"), 5, 100);
          ImGui::SliderInt("maxcount", &param<int>("maxcount"), 100, 2000);

          const char* items[] = { "bisector", "radius","thirdopt" };
          static const char* item_current = items[0];            // Here our selection is a single pointer stored outside the object.
          if (ImGui::BeginCombo("segmentation method", item_current)) // The second parameter is the label previewed before opening the combo.
          {
              for (int n = 0; n < IM_ARRAYSIZE(items); n++)
              {
                  bool is_selected = (item_current == items[n]);
                  if (ImGui::Selectable(items[n], is_selected))
                      item_current = items[n];
                  if (is_selected)
                      ImGui::SetItemDefaultFocus();
                  if (item_current == items[0])
                      current_method = masb::bisector;
                  else if (item_current == items[1])
                      current_method = masb::radius;
                  else if (item_current == items[2])
                      current_method = masb::thirdopt;
                  // Set the initial focus when opening the combo (scrolling + for keyboard navigation support in the upcoming navigation branch)
              }
              ImGui::EndCombo();
          }

          if (item_current == items[0]) {
              //ImGui::SliderFloat("bisec_thres", &params.bisec_thres, 0, 90);
              ImGui::SliderFloat("bisec_thres", &param<float>("bisec_thres"), 0, 90);
          }
          else if (item_current == items[1]) {
              ImGui::SliderFloat("balloverlap_thres", &param<float>("balloverlap_thres"), 0, 100);
              ImGui::Text("it is time for radius");
          }
          else if (item_current == items[2]) {
              ImGui::Text("it is time for third opt");
          }
      }
      void process();
  };
  class MaPt_in_oneTraceNode :public Node {
  public:
      using Node::Node;
      void init() {
          add_param("SearchRadius", (float) 45.00);
          add_input("madata", TT_any);
          add_input("maGeometry", TT_any);
          //add_input("segmentation", TT_any);
          add_input("sheets", TT_any);

          add_output("candidate_r",TT_point_collection);
          add_output("candidate_cos", TT_point_collection);
          add_output("lable",TT_vec1i);
          add_output("traces", TT_line_string_collection);
      }
      void process();
  };
  class ExtractCandidatePtNode :public Node {
  public:
      using Node::Node;
      void init() {
          add_param("SearchRadius", (float) 45.00);
          add_input("madata", TT_any);
          add_input("maGeometry", TT_any);
          //add_input("segmentation", TT_any);
          add_input("sheets", TT_any);

          add_output("candidate_r", TT_point_collection);
          add_output("candidate_cos", TT_point_collection);
      }
      void process();
  };

  class ReadCandidatePtNode :public Node {
  public:
      using Node::Node;
      void init() {
          add_output("candidate_r", TT_point_collection);
          add_output("directon", TT_vec3f);
          add_output("seg_id", TT_vec1i);
      }
      void process();
  };

  class ReadCandidatePtWithBisecNode :public Node {
  public:
      using Node::Node;
      void init() {
          add_output("candidate_r", TT_point_collection);
          add_output("directon", TT_vec3f);
          add_output("seg_id", TT_vec1i);
          add_output("bisector_p", TT_vec3f);
          add_output("bisector_q", TT_vec3f);
      }
      void process();
  };

  class ConnectCandidatePtNode :public Node {
  public:
      //masb::_pram params;
      using Node::Node;
      void init() {
          add_input("pointCloud", TT_point_collection);
          add_input("candidate", TT_point_collection);
          add_input("directon", TT_vec3f);
          add_input("seg_id", TT_vec1i);
          add_input("bisector_p", TT_vec3f);
          add_input("bisector_q", TT_vec3f);
          add_output("bisector_p_vis", TT_line_string_collection);
          add_output("bisector_q_vis", TT_line_string_collection);
          add_output("directon_vis", TT_line_string_collection);
          add_output("directon2_vis", TT_line_string_collection);
          add_output("filter", TT_vec1i);
      }
      void process();
  };
  class PLYLoaderNode :public Node {
  public:
      using Node::Node;
      void init() {
          add_param("filepath", (std::string) "C:/Users/wangq/Downloads/thesis/p3_data/OriginalPointCloud.ply");
          add_param("thinning_factor", (int)50);

          add_output("PointCloud", TT_point_collection);
      }

      void gui() {
          ImGui::InputText("filepath", &param<std::string>("filepath"));
          ImGui::SliderInt("thinning_factor", &param<int>("thinning_factor"), 1, 200);
      }
      void process();
  };

}