#include "ConnectCandidatePt.h"
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/kruskal_min_spanning_tree.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <iostream>
#include <fstream>
#include <map>

using namespace boost;

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
    masb::intList &seg_id, masb::intList &filter, ridge::segment &segmentList, masb::intList &idList,
    line &symple_segmentList, masb::intList &symple_idList){
    ////////////////////////////////////////
    float filter_thresh = 100.0;
    ////////////////////////////////////////
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
        if (cur_sheet == 0)
            continue;
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
        float maxWeight = 0;
        int maxWeight_idx1, maxWeight_idx2;
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
                auto wi = Vrui::Geometry::dist(cur_candidate[i], cur_candidate[j]);
                weightmap[e] = wi;
                if (wi > maxWeight) {
                    maxWeight = wi;
                    maxWeight_idx1 = i;
                    maxWeight_idx2 = j;
                }
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
        
        
        Graph g_short;
        for (std::size_t j = 0; j < cur_size; ++j) {
            Vertex vdi = boost::add_vertex(Vertex_Id(j), g_short);
        }
        property_map<Graph, edge_weight_t>::type weightmap_s = get(edge_weight, g_short);
        for (std::vector < Edge >::iterator ei = spanning_tree.begin();
            ei != spanning_tree.end(); ++ei) {
            auto idx_p = source(*ei, g);
            auto idx_q = target(*ei, g);
            E ea = E(idx_p, idx_q);
            Edge e; bool inserted;
            boost::tie(e, inserted) = add_edge(ea.first, ea.second, g_short);
            weightmap_s[e] = weight[*ei];
        }


        std::vector<Vertex> p(num_vertices(g_short));
        std::vector<float> d(num_vertices(g_short));
        property_map<Graph, vertex_index_t>::type indexmap = get(vertex_index, g_short);
        //property_map<Graph, edge_weight_t>::type weightmap_s = get(edge_weight, g_short);



        /*float maxWeight = 0;
        Vertex maxEnd;
        Vertex maxStart;
        */
        //for (std::size_t j = 0; j < cur_size; ++j) {
            //Vertex s = vertex(j, g_short);
            //dijkstra_shortest_paths(g_short, s, &p[0], &d[0], weightmap_s, indexmap,
            //    std::less<int>(), closed_plus<int>(),
            //    (std::numeric_limits<int>::max)(), 0,
            //    default_dijkstra_visitor());
            //dijkstra_shortest_paths(g_short, s, predecessor_map(&p[0]).distance_map(&d[0]));
            /*
            std::cout << "distances and parents:" << std::endl;
            graph_traits < Graph >::vertex_iterator vi, vend;
            for (tie(vi, vend) = vertices(g_short); vi != vend; ++vi) {
                std::cout << "distance(" << indexmap[*vi] << ") = " << d[*vi] << ", ";
                std::cout << "parent(" << indexmap[*vi] << ") = " << indexmap[p[*vi]] << std::
                    endl;
            }
            std::cout << std::endl;
            */
        /*
            graph_traits < Graph >::vertex_iterator vii, vendd;
            for (tie(vii, vendd) = vertices(g_short); vii != vendd; ++vii) {
                if (maxWeight < d[*vii]) {
                    maxWeight = d[*vii];
                    maxEnd = *vii;
                    maxStart = s;
                }
            }           
        }
        
        
        //because d,p are modified, need to calculate shortest path again
        */

        //dijkstra_shortest_paths(g_short, maxStart, predecessor_map(&p[0]).distance_map(&d[0]));
       

        Vertex maxStart = vertex(maxWeight_idx1, g_short);
        Vertex maxEnd = vertex(maxWeight_idx2, g_short);

        std::cout << "longest path -- " << maxWeight
            << " starts at point -- " << indexmap[maxStart]<<" ("<< maxWeight_idx1<<")"
            << " ends at point -- " << indexmap[maxEnd] << " (" << maxWeight_idx2 << ")\n";

        dijkstra_shortest_paths(g_short, maxStart, predecessor_map(&p[0]).distance_map(&d[0]));

        std::vector< graph_traits< Graph >::vertex_descriptor > path;
        graph_traits< Graph >::vertex_descriptor current = maxEnd;
        
        while (current != maxStart) {
            path.push_back(current);
            current = p[current];
        }
        path.push_back(maxStart);
        masb::PointList a_line;
        std::vector< graph_traits< Graph >::vertex_descriptor >::iterator it;
        for (it = path.begin(); it != path.end(); ++it) {
            //std::cout << indexmap[*it] << " ";
            auto pt_idx = indexmap[*it];
            masb::Point pt = cur_candidate[pt_idx];
            //std::cout << pt[0] << " " << pt[1] << " " << pt[2] << std::endl;
            a_line.push_back(pt);
        }
        //std::cout << std::endl;
        symple_segmentList.push_back(a_line);
        symple_idList.push_back(cur_sheet);



        /*
        //example for shortest path
        typedef adjacency_list < listS, vecS, directedS,no_property, property < edge_weight_t, float > > graph_t;
        typedef graph_traits < graph_t >::vertex_descriptor vertex_descriptor;
        typedef graph_traits < graph_t >::edge_descriptor edge_descriptor;
        typedef std::pair<int, int> Edge;

        const int num_nodes = 5;
        enum nodes { A, B, C, D, E };
        char name[] = "ABCDE";
        Edge edge_array[] = { Edge(A, C), Edge(B, B), Edge(B, D), Edge(B, E),
          Edge(C, B), Edge(C, D), Edge(D, E), Edge(E, A), Edge(E, B)
        };
        float weights[] = { 1.1, 2.1, 1.2, 2.2, 7.7, 3.4, 1.6, 1.76, 1.0 };
        int num_arcs = sizeof(edge_array) / sizeof(Edge);
        #if defined(BOOST_MSVC) && BOOST_MSVC <= 1300
        graph_t g(num_nodes);
        property_map<graph_t, edge_weight_t>::type weightmap = get(edge_weight, g);
        for (std::size_t j = 0; j < num_arcs; ++j) {
            edge_descriptor e; bool inserted;
            tie(e, inserted) = add_edge(edge_array[j].first, edge_array[j].second, g);
            weightmap[e] = weights[j];
        }
        #else
        graph_t g(edge_array, edge_array + num_arcs, weights, num_nodes);
        property_map<graph_t, edge_weight_t>::type weightmap = get(edge_weight, g);
        #endif
        std::vector<vertex_descriptor> p(num_vertices(g));
        std::vector<float> d(num_vertices(g));
        vertex_descriptor s = vertex(A, g);

        #if defined(BOOST_MSVC) && BOOST_MSVC <= 1300
        // VC++ has trouble with the named parameters mechanism
        property_map<graph_t, vertex_index_t>::type indexmap = get(vertex_index, g);
        dijkstra_shortest_paths(g, s, &p[0], &d[0], weightmap, indexmap,
            std::less<int>(), closed_plus<int>(),
            (std::numeric_limits<int>::max)(), 0,
            default_dijkstra_visitor());
        #else
        dijkstra_shortest_paths(g, s, predecessor_map(&p[0]).distance_map(&d[0]));
        #endif

        std::cout << "distances and parents:" << std::endl;
        graph_traits < graph_t >::vertex_iterator vi, vend;
        for (tie(vi, vend) = vertices(g); vi != vend; ++vi) {
            std::cout << "distance(" << name[*vi] << ") = " << d[*vi] << ", ";
            std::cout << "parent(" << name[*vi] << ") = " << name[p[*vi]] << std::
                endl;
        }
        std::cout << std::endl;

        */
        
    }
}
void ridge::connectCandidatePt8MST_nosegid(masb::PointList &PointCloud, masb::PointList & candidate, masb::intList &filter, segment &segmentList) {
    /////////////////////////////////////////
    float filter_thresh = 100.0;
    float isolate_thresh = 50.0;
    float length_thresh = 100.0;
    /////////////////////////////////////////
    std::cout << "before filter candidate pt -- " << candidate.size() << "\n";
    kdtree2::KDTree* pc_kdtree;
    pc_kdtree = new kdtree2::KDTree(PointCloud, true);
    pc_kdtree->sort_results = true;

    kdtree2::KDTree* candidate_kdtree;
    candidate_kdtree = new kdtree2::KDTree(candidate, true);
    candidate_kdtree->sort_results = true;

    std::vector<float> dis_list; dis_list.reserve(candidate.size());
    for (int i = candidate.size() - 1; i >= 0; i--) {
        kdtree2::KDTreeResultVector neighbours, neighbours_can;
        pc_kdtree->n_nearest(candidate[i], 1, neighbours);
        dis_list.push_back(neighbours[0].dis);

        candidate_kdtree->n_nearest(candidate[i], 2, neighbours_can);

        if (neighbours[0].dis > filter_thresh ||
            neighbours_can[1].dis>isolate_thresh) {
            filter.push_back(0);
            candidate.erase(candidate.begin() + i);
        }
        else
            filter.push_back(1);
    }

    std::cout << "after filter candidate pt -- " << candidate.size() << "\n";

    size_t cur_size = candidate.size();
    Graph g(cur_size);
    size_t num_edges = cur_size * (cur_size - 1) / 2;
    property_map<Graph, edge_weight_t>::type weightmap = get(edge_weight, g);
    //property_map<Graph, vertex_index_t>::type VertexIndexMap = get(vertex_index, g);
    /////////////////////////////////////////////////////////////
    //       large dataset have unknown bug here!!!!!!!!!!!!!!
    /////////////////////////////////////////////////////////////
    for (size_t i = 0; i < cur_size; ++i) {
        for (size_t j = i + 1; j < cur_size; ++j) {
            E ea = E(i, j);
            Edge e; bool inserted;
            boost::tie(e, inserted) = add_edge(ea.first, ea.second, g);
            auto wi = Vrui::Geometry::sqrDist(candidate[i], candidate[j]);
            weightmap[e] = wi;
            std::cout << "i,j---" << i << "--" << j << std::endl;
        }
    }

    property_map < Graph, edge_weight_t >::type weight = get(edge_weight, g);
    std::vector < Edge > spanning_tree;

    auto tmp = std::back_inserter(spanning_tree);
    kruskal_minimum_spanning_tree(g, tmp);
    for (std::vector < Edge >::iterator ei = spanning_tree.begin();
        ei != spanning_tree.end(); ++ei) {
        if (weight[*ei] < length_thresh) {
            auto idx_p = source(*ei, g);
            auto idx_q = target(*ei, g);
            masb::Point p = candidate[idx_p];
            masb::Point q = candidate[idx_q];
            segmentList.push_back(ridge::PointPair(p, q));
        }
    }
}