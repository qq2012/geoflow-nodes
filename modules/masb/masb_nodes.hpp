#include <geoflow/core/geoflow.hpp>
#include <geoflow/gui/osdialog.hpp>

#include <compute_ma_processing.h>
#include <compute_normals_processing.h>
#include "Ma_geometry_processing.hpp"
#include "Ma_utility.h"
#include "Ma_Segmentation_processing.h"
#include "sheet_adjacency.h"
#include "ExtractCandidatePt.h"
#include "MaPt_in_oneTrace.h"
#include "ConnectCandidatePt.h"
#include "breaklineValidate.h"


namespace geoflow::nodes::mat {

    class ComputeMedialAxisNode :public Node {
    public:
        float interval = 2;
        double zero = 0, pi = 3.14;
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
        void gui() {
            ImGui::SliderFloat("initial_radius", &param<float>("initial_radius"), 0, 1000);
            ImGui::SliderScalar("denoise_preserve", ImGuiDataType_Double, &param<double>("denoise_preserve"), &zero, &pi);
            ImGui::SliderScalar("denoise_planar", ImGuiDataType_Double, &param<double>("denoise_planar"), &zero, &pi);
            ImGui::Checkbox("nan_for_initr", &param<bool>("nan_for_initr"));
        }
        void process();
    };

    class ComputeNormalsNode :public Node {
    public:
        float interval = 2;
        using Node::Node;
        void init() {
            add_param("K", (int)10);
            add_input("points", typeid(PointCollection));
            add_output("normals", typeid(vec3f));
        }
        void gui() {
            ImGui::SliderInt("K", &param<int>("K"), 1, 100);
        }
        void process();

    };

    class MaGeometryNode :public Node {
    public:
        using Node::Node;
        void init() {
            add_param("isInterior", (bool)false);// 

            add_input("points", typeid(PointCollection));
            add_input("normals", typeid(vec3f));
            add_input("ma_coords", typeid(PointCollection));
            add_input("ma_qidx", typeid(vec1i));
            add_input("ma_radius", typeid(vec1f));
            
            add_output("atom", typeid(PointCollection));
            add_output("unshrikingGroundPoint", typeid(PointCollection));
            add_output("sp_reverse_norm", typeid(vec3f));
            add_output("spokeVectorP", typeid(vec3f));
            //add_output("spokeVectorQ", typeid(vec3f));
            add_output("bisector", typeid(vec3f));
            add_output("ma_direction", typeid(vec3f));
            add_output("mat", typeid(masb::MAT));
        }
        void gui() {
            ImGui::Checkbox("isInterior", &param<bool>("isInterior"));
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
            add_param("mincount", (int)50);
            add_param("maxcount", (int)1000);
            add_param("bisec_thres", (float)3.0);
            add_param("spokecross_thres", (float)3.0);
            add_param("balloverlap_thres", (float)5.0);

            add_input("mat", typeid(masb::MAT));
            add_output("seg_id", typeid(vec1i));
            add_output("sheets", typeid(masb::Sheet_idx_List));
        }
        void gui() {
            ImGui::SliderInt("mincount", &param<int>("mincount"), 10, 200);
            ImGui::SliderInt("maxcount", &param<int>("maxcount"), 100, 2000);

            const char* items[] = { "bisector", "spokecross","balloverlap",
                "combinBisecAndSpcros","combinBallAndSpcros" };

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
                        current_method = masb::spokecross;
                    else if (item_current == items[2])
                        current_method = masb::balloverlap;
                    else if (item_current == items[3])
                        current_method = masb::combinBisecAndSpcros;
                    else if (item_current == items[4])
                        current_method = masb::combinBallAndSpcros;
                    // Set the initial focus when opening the combo (scrolling + for keyboard navigation support in the upcoming navigation branch)
                }
                ImGui::EndCombo();
            }

            if (item_current == items[0]) {
                ImGui::SliderFloat("bisec_thres", &param<float>("bisec_thres"), 0, 45);
            }
            else if (item_current == items[1]) {
                ImGui::SliderFloat("spokecross_thres", &param<float>("spokecross_thres"), 0, 45);
            }
            else if (item_current == items[2]) {
                ImGui::SliderFloat("balloverlap_thres", &param<float>("balloverlap_thres"), 0, 45);
            }
            else if (item_current == items[3]) {
                ImGui::SliderFloat("bisec_thres", &param<float>("bisec_thres"), 0, 45);
                ImGui::SliderFloat("spokecross_thres", &param<float>("spokecross_thres"), 0, 45);
            }
            else if (item_current == items[4]) {
                ImGui::SliderFloat("balloverlap_thres", &param<float>("balloverlap_thres"), 0, 10);
                ImGui::SliderFloat("spokecross_thres", &param<float>("spokecross_thres"), 0, 45);
            }
        }
        void process();
    };
    class adjacencyNode :public Node {
    public:
        using Node::Node;
        void init() {
            add_param("searchRadius", (float) 20.00);
            add_param("adjacency_thresh", (int)10);
            add_input("mat", typeid(masb::MAT));
            add_input("sheets", typeid(masb::Sheet_idx_List));
            add_output("sheet-sheet adjacency", typeid(ridge::int_pair_vec));
            add_output("junction points", typeid(vec3f));
        }
        void gui() {
            ImGui::SliderFloat("searchRadius", &param<float>("searchRadius"), 0, 200);
            ImGui::SliderInt("adjacency_thresh", &param<int>("adjacency_thresh"), 0, 100);
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

            add_output("candidate_r", typeid(PointCollection));
            add_output("candidate_cos", typeid(PointCollection));
            add_output("lable", typeid(vec1i));
            add_output("traces", typeid(LineStringCollection));
        }

        void process();
    };
    class ExtractCandidatePtNode :public Node {
    public:
        using Node::Node;
        void init() {
            add_param("SearchRadius", (float) 45.00);
            add_param("deviationAng_thres", (float) 15.00);
            add_param("filterDistance", (float) 5.00);
            add_param("bis_avg_knn", (int)20);

            add_input("mat", typeid(masb::MAT));
            add_input("seg_id", typeid(vec1i));
            add_input("pointcloud", typeid(PointCollection));

            //add_output("edgeBall", typeid(masb::MAT));
            add_output("edgeBallAtom", typeid(PointCollection));
            //add_output("sp_reverse_norm", typeid(vec3f));
            //add_output("spokeVectorP", typeid(vec3f));
            //add_output("spokeVectorQ", typeid(vec3f));
            //add_output("bisector", typeid(vec3f));
            //add_output("ma_direction", typeid(vec3f));
            add_output("edgeBall_id", typeid(vec1i));

            add_output("candidate_r", typeid(PointCollection));
            add_output("candidate_cos", typeid(PointCollection));
            add_output("candidate_r_bisector_avg", typeid(PointCollection));

            add_output("filter", typeid(vec1i));
            add_output("candidate_points", typeid(PointCollection));
            add_output("candidate_points_id", typeid(vec1i));
            add_output("candidate_points_radius", typeid(vec1f));
            add_output("candidate_points_direction", typeid(vec3f));
        }
        void gui() {
            ImGui::SliderFloat("SearchRadius", &param<float>("SearchRadius"), 20, 100);
            ImGui::SliderFloat("deviationAng_thres", &param<float>("deviationAng_thres"), 0, 45);
            ImGui::SliderInt("bis_avg_knn", &param<int>("bis_avg_knn"), 5, 50);
            ImGui::SliderFloat("filterDistance", &param<float>("filterDistance"), 0.1, 10);
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
            add_param("filepath", (std::string) "C:/Users/wangq/Downloads/thesis/p3_data/candidatept_smalltest2_r100SepAng15np500_Bis_avg_all_WithBisector.ply");

            add_output("candidate_r", typeid(PointCollection));
            add_output("directon", typeid(vec3f));
            add_output("seg_id", typeid(vec1i));
            add_output("bisector_p", typeid(vec3f));
            add_output("bisector_q", typeid(vec3f));
        }
        void gui() {
            ImGui::InputText("filepath", &param<std::string>("filepath"));
        }
        void process();
    };
    class ReadSegmentRestltNode :public Node {
    public:

        using Node::Node;
        void init() {
            add_param("filepath", (std::string) "C:/Users/wangq/Downloads/thesis/p3_data/samllTest2_segmentResult.ply");

            add_output("ma_coords", typeid(PointCollection));
            add_output("seg_id", typeid(vec1i));
            add_output("radius", typeid(vec1f));
            add_output("cp", typeid(vec3f));
            add_output("cq", typeid(vec3f));
            add_output("bisector", typeid(vec3f));
            add_output("directon", typeid(vec3f));
            add_output("theta", typeid(vec1f));
            add_output("direction_vis", typeid(LineStringCollection));
        }
        void gui() {
            ImGui::InputText("filepath", &param<std::string>("filepath"));
        }
        void process();
    };
    class ReadJunctionPtNode :public Node {
    public:
        using Node::Node;
        void init() {
            add_param("filepath", (std::string) "C:/Users/wangq/Downloads/thesis/p3_data/junctionPt_smalltest2_.ply");

            add_output("ma_coords", typeid(PointCollection));
            add_output("sheet-sheet", typeid(ridge::int_pair_vec));
        }
        void gui() {
            ImGui::InputText("filepath", &param<std::string>("filepath"));
        }
        void process();
    };
    class ReadAdjacencyNode : public Node {
    public:
        using Node::Node;
        void init() {
            add_param("filepath", (std::string) "C:/Users/wangq/Downloads/thesis/p3_data/Adjacency_smalltest2_.txt");

            add_output("sheet-sheet", typeid(ridge::int_pair_vec));
        }
        void gui() {
            ImGui::InputText("filepath", &param<std::string>("filepath"));
        }
        void process();

    };

    class SegmentRestlt_TraceNode :public Node {
    public:

        using Node::Node;
        void init() {
            add_input("ma_coords", typeid(PointCollection));
            add_input("seg_id", typeid(vec1i));
            add_input("radius", typeid(vec1f));
            add_input("cp", typeid(vec3f));
            add_input("cq", typeid(vec3f));
            add_input("bisector", typeid(vec3f));
            add_input("directon", typeid(vec3f));
            add_input("theta", typeid(vec1f));

            add_output("madata", typeid(masb::mat_data));
            add_output("maGeometry", typeid(masb::ma_Geometry));
            add_output("sheets", typeid(masb::Sheet_idx_List));
        }
        void process();
    };
    class ConnectCandidatePtNode :public Node {
    public:
        using Node::Node;
        void init() {
            add_input("candidate_points", typeid(PointCollection));
            add_input("point_directon", typeid(vec3f));
            add_input("candidate_points_id", typeid(vec1i));
            add_input("adjacency", typeid(ridge::int_pair_vec));

            add_output("mstLineSegment", typeid(LineStringCollection));
            add_output("mstLineSegment_id", typeid(vec1i));
            add_output("polylines_maxDistance", typeid(LineStringCollection));
            add_output("polylines_maxAccDist", typeid(LineStringCollection));
            //output("polylines_maxPts", typeid(LineStringCollection));
            add_output("polyline_id", typeid(vec1i));

            add_output("linesWithJunction_maxAccDist", typeid(LineStringCollection));
            add_output("linesWithJunction_maxDistance", typeid(LineStringCollection));
        }
        void process();
    };
    class PLYLoaderNode :public Node {
    public:
        using Node::Node;
        void init() {
            //add_param("filepath", (std::string) "C:/Users/wangq/Downloads/thesis/p3_data/OriginalPointCloud.ply");
            add_param("filepath", (std::string) "C:/Users/wangq/Downloads/thesis/p3_data/smallTest2_pointCloud.ply");
            add_param("thinning_factor", (int)50);

            add_output("PointCloud", typeid(PointCollection));
        }

        void gui() {
            ImGui::InputText("filepath", &param<std::string>("filepath"));
            ImGui::SliderInt("thinning_factor", &param<int>("thinning_factor"), 1, 200);
        }
        void process();
    };
    class LoadReferenceBreaklineNode :public Node {
    public:
        using Node::Node;
        void init() {
            add_output("ReferenceVertices", typeid(PointCollection));
            add_output("ReferenceBreakline", typeid(LineString));

            add_param("filepath", (std::string) "");
        }
        void gui() {
            ImGui::FilePicker(OSDIALOG_OPEN, param<std::string>("filepath"));
        }
        void process();
    };

    class LoadTruePositiveVerticesNode : public Node {
    public:
        using Node::Node;
        void init() {
            add_output("ExtractedTruePositiveVertices", typeid(PointCollection));
            add_output("ExtractedTruePositiveBreakline", typeid(LineString));

            add_param("filepath", (std::string) "");
        }
        void gui() {
            ImGui::FilePicker(OSDIALOG_OPEN, param<std::string>("filepath"));
        }
        void process();

    };

    class SelectTestBreaklineNode :public Node {
    public:
        using Node::Node;
        void init() {
            add_param("BreaklineID", (int)0);
            add_input("ExtractedBreakline", typeid(LineStringCollection));

            add_output("TestBreakline", typeid(LineString));
        }
        void gui() {
            ImGui::SliderInt("BreaklineID", &param<int>("BreaklineID"), 1, 20);
        }
        void process();
    };

    class BreaklineValidationNode :public Node {
    public:
        using Node::Node;
        void init() {
            add_input("ReferenceVertices", typeid(PointCollection));
            add_input("TruePositiveVertices", typeid(PointCollection));//"TestBreakline"
        }
        void process();
    };

}