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
        add_input("points", typeid(PointCollection));
        add_input("normals", typeid(vec3f));
        add_output("ma_coords", typeid(PointCollection));
        add_output("ma_radii", typeid(vec1f));
        add_output("ma_radius", typeid(vec1f));
        add_output("ma_qidx", typeid(vec1i));
        add_output("ma_is_interior", typeid(vec1i));
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
    float interval = 2;
    using Node::Node;
    void init() {
      add_param("K", (int) 10);
      add_input("points", typeid(PointCollection));
      add_output("normals", typeid(vec3f));
    }
    void gui(){
      ImGui::SliderInt("K", &param<int>("K"), 1, 100);
    }
    void process();

  };

  class MaGeometryNode :public Node {
  public:
      using Node::Node;
      void init() {
          add_input("points", typeid(PointCollection));
          add_input("normals", typeid(vec3f));
          add_input("ma_coords", typeid(PointCollection));
          add_input("ma_qidx", typeid(vec1i));

          add_output("ma_SeparationAng", typeid(vec1f));
          add_output("ma_bisector", typeid(vec3f));
          add_output("bisec_vis",typeid(LineStringCollection));
          add_output("ma_normal", typeid(vec3f));
          add_output("ma_normal_vis", typeid(LineStringCollection));
      }
      void process();
  };
  class FilterRNode :public Node {
  public:
      //float radius = 200.00;
      using Node::Node;
      void init() {
          add_param("filterRadius", (float) 100.00);
          add_input("ma_radius", typeid(vec1f));
          add_output("remaining_idx", typeid(vec1i));
      }
      void gui() {
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

          add_input("remaining_idx", typeid(vec1i));
          add_input("ma_coords", typeid(PointCollection));
          add_input("ma_qidx", typeid(vec1i));
          add_input("ma_radius", typeid(vec1f));
          add_input("ma_SeparationAng", typeid(vec1f));
          add_input("ma_bisector", typeid(vec3f));
          add_input("ma_normal", typeid(vec3f));

          add_output("seg_id", typeid(vec1i));
          add_output("sheet_all", typeid(masb::Sheet_idx_List));//#########  todo type -- line string clooection
          add_output("ma_coords", typeid(PointCollection));
          add_output("madata_in", typeid(masb::mat_data));
          add_output("madata_out", typeid(masb::mat_data));
          add_output("maGeometry_in",typeid(masb::ma_Geometry));
          add_output("maGeometry_out", typeid(masb::ma_Geometry));
          add_output("seg_in", typeid(std::vector<long long int>));
          add_output("seg_out", typeid(std::vector<long long int>));
          add_output("sheet_in", typeid(masb::Sheet_idx_List));
          add_output("sheet_out", typeid(masb::Sheet_idx_List));
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
          add_input("madata", typeid(masb::mat_data));
          add_input("maGeometry", typeid(masb::ma_Geometry));
          add_input("sheets", typeid(masb::Sheet_idx_List));

          add_output("candidate_r",typeid(PointCollection));
          add_output("candidate_cos", typeid(PointCollection));
          add_output("lable",typeid(vec1i));
          add_output("traces", typeid(LineStringCollection));
      }
      void process();
  };
  class ExtractCandidatePtNode :public Node {
  public:
      using Node::Node;
      void init() {
          add_param("SearchRadius", (float) 45.00);
          add_input("madata", typeid(masb::mat_data));
          add_input("maGeometry", typeid(masb::ma_Geometry));
          add_input("sheets", typeid(masb::Sheet_idx_List));

          add_output("candidate_r", typeid(PointCollection));
          add_output("candidate_cos", typeid(PointCollection));
      }
      void process();
  };

  class ReadCandidatePtNode :public Node {
  public:
      using Node::Node;
      void init() {
          add_output("candidate_r", typeid(PointCollection));
          add_output("directon", typeid(vec3f));
          add_output("seg_id", typeid(vec1i));
      }
      void process();
  };

  class ReadCandidatePtWithBisecNode :public Node {
  public:
      using Node::Node;
      void init() {
          add_output("candidate_r", typeid(PointCollection));
          add_output("directon", typeid(vec3f));
          add_output("seg_id", typeid(vec1i));
          add_output("bisector_p", typeid(vec3f));
          add_output("bisector_q", typeid(vec3f));
      }
      void process();
  };

  class ConnectCandidatePtNode :public Node {
  public:
      //masb::_pram params;
      using Node::Node;
      void init() {
          add_input("pointCloud", typeid(PointCollection));
          add_input("candidate", typeid(PointCollection));
          add_input("directon", typeid(vec3f));
          add_input("seg_id", typeid(vec1i));
          add_input("bisector_p", typeid(vec3f));
          add_input("bisector_q", typeid(vec3f));

          add_output("filter", typeid(vec1i));
          add_output("ridge", typeid(LineStringCollection));
          add_output("ridgeId", typeid(vec1i));
          add_output("longest_path", typeid(LineStringCollection));
          add_output("longest_id",typeid(vec1i));
          add_output("smoothLine", typeid(LineStringCollection));
          add_output("smoothpoint", typeid(PointCollection));

          add_output("bisector_p_vis", typeid(LineStringCollection));
          add_output("bisector_q_vis", typeid(LineStringCollection));
          add_output("directon_vis", typeid(LineStringCollection));
          add_output("directon2_vis", typeid(LineStringCollection));
          
      }
      void process();
  };
  class PLYLoaderNode :public Node {
  public:
      using Node::Node;
      void init() {
          add_param("filepath", (std::string) "C:/Users/wangq/Downloads/thesis/p3_data/OriginalPointCloud.ply");
          add_param("thinning_factor", (int)50);

          add_output("PointCloud", typeid(PointCollection));
      }

      void gui() {
          ImGui::InputText("filepath", &param<std::string>("filepath"));
          ImGui::SliderInt("thinning_factor", &param<int>("thinning_factor"), 1, 200);
      }
      void process();
  };

}