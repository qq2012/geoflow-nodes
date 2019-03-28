#include "ConnectCandidatePt.h"
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/kruskal_min_spanning_tree.hpp>
#include <iostream>
#include <fstream>

struct VertexId { int idx; };
using namespace boost;
//typedef adjacency_list < vecS, vecS, undirectedS, no_property, 
//                                property < edge_weight_t, float > > Graph;
typedef property<boost::edge_weight_t, float> EdgeWeightProperty;
typedef property<boost::vertex_index_t, size_t> Vertex_Id;
typedef adjacency_list<vecS, vecS, undirectedS, Vertex_Id, EdgeWeightProperty> Graph;
typedef graph_traits < Graph >::edge_descriptor Edge;
typedef graph_traits < Graph >::vertex_descriptor Vertex;//VertexDescriptor
typedef std::pair<int, int> E;
typedef boost::graph_traits<Graph>::vertex_iterator     VertexIterator;
// edge iterator (from which you can traverse all the edges of a graph)
//typedef boost::graph_traits<Graph>::edge_iterator       EdgeIterator;

void ridge::connectCandidatePt8MST(masb::PointList &PointCloud, masb::PointList &candidate,
    masb::intList &seg_id, masb::intList &filter, ridge::segment &segmentList, masb::intList &idList){
    
    float filter_thresh = 5.0;

    kdtree2::KDTree* pc_kdtree;
    pc_kdtree = new kdtree2::KDTree(PointCloud, true);
    pc_kdtree->sort_results = true;
    std::vector<float> dis_list; dis_list.reserve(candidate.size());
    for (int i = candidate.size() - 1; i >= 0; i--) {
        kdtree2::KDTreeResultVector neighbours;
        pc_kdtree->n_nearest(candidate[i], 1, neighbours);
        dis_list.push_back(neighbours[0].dis);
        if (neighbours[0].dis > filter_thresh) {
            filter.push_back(0);
            candidate.erase(candidate.begin() + i);
            seg_id.erase(seg_id.begin() + i);
        }
        else
            filter.push_back(1);
    }
    //dis_list;

    std::map<int, int> seg_frequency;
    for (int i : seg_id)
        ++seg_frequency[i];
    //for (const auto& e : seg_frequency)
    //    std::cout << "Element " << e.first 
    //     << " encountered " << e.second << " times\n";

    for (const auto& e : seg_frequency) {
        auto cur_sheet = e.first;
        auto cur_size = e.second;
        masb::PointList cur_candidate;
        cur_candidate.reserve(cur_size);
        for (int i = 0; i < seg_id.size(); i++) {
            if (seg_id[i] == cur_sheet) {
                cur_candidate.push_back(candidate[i]);
            }
        }

        if (cur_candidate.size() != cur_size)
            std::cout << "Error" << std::endl;

        /*
        //example
        const int num_nodes = 5;
        E edge_array[] = { E(0, 2), E(1, 3), E(1, 4), E(2, 1), E(2, 3),
          E(3, 4), E(4, 0), E(4, 1)
        };
        int weights[] = { 1, 1, 2, 7, 3, 1, 1, 1 };
        std::size_t num_edges = sizeof(edge_array) / sizeof(E);
        Graph g(edge_array, edge_array + num_edges, weights, num_nodes);
        */
        /*
        //initialize graph  -- wrong??why????????????????????????
        cur_size = 20;
        Graph g;
        for (int i = 0; i < cur_size; i++) {
            Vertex vdi = boost::add_vertex(Vertex_Id(i), g);
        }
        if (boost::num_vertices(g) != cur_size)
            std::cout << "error in adding vertices" << std::endl;

        std::pair<VertexIterator, VertexIterator> v = boost::vertices(g);
        int i = 0;
        for (VertexIterator vi = v.first; vi != v.second; ++vi) {
            //vi.first -- begin ; vi.second -- end;
            Vertex vdi = *vi;
            int j = i + 1;
            for (VertexIterator vj = vi+1; vj != v.second; ++vj) {
                Vertex vdj = *vj;
                auto wij = Vrui::Geometry::sqrDist(cur_candidate[i], cur_candidate[j]);
                boost::add_edge(vdi, vdi, EdgeWeightProperty(wij), g);
                //std::cout << "add edge " << i << " -- " << j << std::endl;
                j++;
            }
            i++;
        }
        if (boost::num_edges(g) != cur_size * (cur_size - 1) / 2)
            std::cout << "error in adding edge" << std::endl;
        */
        std::cout << "start connect ridge " << cur_sheet << std::endl;
        //cur_size = 200;
        Graph g(cur_size);
        size_t num_edges = cur_size * (cur_size - 1) / 2;
        property_map<Graph, edge_weight_t>::type weightmap = get(edge_weight, g);
        //property_map<Graph, vertex_index_t>::type VertexIndexMap = get(vertex_index, g);
        for (size_t i = 0; i < cur_size; ++i) {
            for (size_t j = i + 1; j < cur_size; ++j) {
                E ea = E(i, j);
                Edge e; bool inserted;
                boost::tie(e, inserted) = add_edge(ea.first, ea.second, g);
                weightmap[e] =  Vrui::Geometry::sqrDist(cur_candidate[i], cur_candidate[j]);
            }
        }
        
        property_map < Graph, edge_weight_t >::type weight = get(edge_weight, g);
        //property_map < Graph, vertex_index_t >::type index = get(vertex_index, g);
        std::vector < Edge > spanning_tree;

        auto tmp = std::back_inserter(spanning_tree);
        kruskal_minimum_spanning_tree(g, tmp);
        /*
        std::cout << "Print the edges in the MST:" << std::endl;
        for (std::vector < Edge >::iterator ei = spanning_tree.begin();
            ei != spanning_tree.end(); ++ei) {
            std::cout << source(*ei, g) << " <--> " << target(*ei, g)
                << " with weight of " << weight[*ei]
                << std::endl;
        }
        Print the edges in the MST:
        0 <--> 2 with weight of 1
        3 <--> 4 with weight of 1
        4 <--> 0 with weight of 1
        1 <--> 3 with weight of 1
        */
        for (std::vector < Edge >::iterator ei = spanning_tree.begin();
            ei != spanning_tree.end(); ++ei) {
            auto idx_p = source(*ei, g);
            auto idx_q = target(*ei, g);
            masb::Point p = cur_candidate[idx_p];
            masb::Point q = cur_candidate[idx_q];
            segmentList.push_back(ridge::PointPair(p, q));
        }
        masb::intList temp;
        temp.resize(cur_size - 1, cur_sheet);
        idList.insert(idList.end(), temp.begin(), temp.end());
    }
}