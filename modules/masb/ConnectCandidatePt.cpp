#include "ConnectCandidatePt.h"
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/kruskal_min_spanning_tree.hpp>
#include <iostream>
#include <fstream>

struct VertexId { int idx; };


masb::intList ridge::connectCandidatePtProcess(masb::PointList &pointCloud,masb::PointList &candidate, 
    masb::intList &seg_id,masb::VectorList &direction, 
    masb::VectorList &bisec_p, masb::VectorList &bisec_q) {

    masb::intList filter;

    kdtree2::KDTree* pc_kdtree;
    pc_kdtree = new kdtree2::KDTree(pointCloud, true);
    pc_kdtree->sort_results = true;
    std::vector<float> dis_list; dis_list.reserve(candidate.size());
    for (int i = candidate.size() - 1; i >= 0; i--) {
        kdtree2::KDTreeResultVector neighbours;
        pc_kdtree->n_nearest(candidate[i], 1, neighbours);
        dis_list.push_back(neighbours[0].dis);
        if (neighbours[0].dis > 30.0) {
            filter.push_back(3);
            /*
            candidate.erase(candidate.begin() + i);
            seg_id.erase(seg_id.begin() + 1);
            direction.erase(direction.begin() + i);
            bisec_p.erase(bisec_p.begin() + i);
            bisec_q.erase(bisec_q.begin() + i);
            */
        }
        else
            filter.push_back(0);
    }
    //dis_list;

    std::map<int, int> seg_frequency;
    for (int i : seg_id)
        ++seg_frequency[i];
    //for (const auto& e : seg_frequency)
    //    std::cout << "Element " << e.first 
    //     << " encountered " << e.second << " times\n";
    for (const auto& e : seg_frequency) {
        auto cur_id = e.first;
        auto cur_size = e.second;
        masb::PointList cur_candidate;
        cur_candidate.reserve(cur_size);
        for (int i = 0; i < seg_id.size(); i++) {
            if (seg_id[i] == cur_id) {
                cur_candidate.push_back(candidate[i]);
            }
        }

        using namespace boost;
        //typedef adjacency_list < vecS, vecS, undirectedS, no_property, 
        //                                property < edge_weight_t, float > > Graph;
        typedef property<boost::edge_weight_t, float> EdgeWeightProperty;
        typedef adjacency_list<vecS, vecS, undirectedS, VertexId, EdgeWeightProperty> Graph;
        typedef graph_traits < Graph >::edge_descriptor Edge;
        typedef graph_traits < Graph >::vertex_descriptor Vertex;//VertexDescriptor
        typedef std::pair<int, int> E;

        int num_nodes = 3;
        E edge_array[] = { E(0, 1), E(0, 2), E(1, 2) };

        auto w1 = Vrui::Geometry::sqrDist(cur_candidate[0], cur_candidate[1]);
        auto w2 = Vrui::Geometry::sqrDist(cur_candidate[0], cur_candidate[2]);
        auto w3 = Vrui::Geometry::sqrDist(cur_candidate[1], cur_candidate[2]);
        float weights[] = { w1,w2,w3 };

        std::size_t num_edges = sizeof(edge_array) / sizeof(E);
        //std::size_t num_edges = num_nodes * (num_nodes - 1) / 2;

#if defined(BOOST_MSVC) && BOOST_MSVC <= 1300
        Graph g(num_nodes);
        property_map<Graph, edge_weight_t>::type weightmap = get(edge_weight, g);
        for (std::size_t j = 0; j < num_edges; ++j) {
            Edge e; bool inserted;
            boost::tie(e, inserted) = add_edge(edge_array[j].first, edge_array[j].second, g);
            weightmap[e] = weights[j];
        }
#else
        Graph g(edge_array, edge_array + num_edges, weights, num_nodes);
#endif
        
        //modify graph
        for (int i = 0; i < num_nodes; i++) {
            for (int j = i + 1; j < num_nodes; j++) {
                if ((i == 0 && j == 1) || (i == 0 && j == 2) || (i == 1 && j == 2))
                    continue;
                VertexId vi, vj;
                vi.idx = i;
                vj.idx = j;
                Vertex vdi = boost::add_vertex(vi, g);
                Vertex vdj = boost::add_vertex(vj, g);
                auto wij = Vrui::Geometry::sqrDist(cur_candidate[i], cur_candidate[j]);
                boost::add_edge(vdi, vdi, EdgeWeightProperty(wij), g);
            }
        };

        property_map < Graph, edge_weight_t >::type weight = get(edge_weight, g);
        std::vector < Edge > spanning_tree;

        kruskal_minimum_spanning_tree(g, std::back_inserter(spanning_tree));

        std::cout << "Print the edges in the MST:" << std::endl;
        for (std::vector < Edge >::iterator ei = spanning_tree.begin();
            ei != spanning_tree.end(); ++ei) {
            std::cout << source(*ei, g) << " <--> " << target(*ei, g)
                << " with weight of " << weight[*ei]
                << std::endl;
        }
    }

    return filter;


}