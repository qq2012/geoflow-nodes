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

void ridge::connectCandidatePt8MST(masb::PointList &pts, masb::VectorList &pt_directon, masb::intList &pt_id,
    ridge::LineSegmentList &mstLineSegment, masb::intList &mstLineSegment_id, 
    ridge::PolyineList &polylines_maxDistance, ridge::PolyineList &polylines_maxAccDist,ridge::PolyineList &polylines_maxPts,
    masb::intList &polyline_id){

    std::map<int, int> seg_frequency;
    for (auto i : pt_id) {
        ++seg_frequency[i];
    }
    for (const auto& e : seg_frequency)
        std::cout << "Element " << e.first 
         << " encountered " << e.second << " times\n";
    for (const auto& e : seg_frequency) {
        auto cur_sheet = e.first;
        auto cur_size = e.second;
        if (cur_sheet == 0)
            continue;
        masb::PointList cur_candidate;
        cur_candidate.reserve(cur_size);
        for (int i = 0; i < pt_id.size(); i++) {
            if (pt_id[i] == cur_sheet) {
                cur_candidate.push_back(pts[i]);
            }
        }
        if (cur_candidate.size() != cur_size)
            std::cout << "Error" << std::endl;
        
        float maxDistance = 0;
        int maxDistance_idx1, maxDistance_idx2;

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
                if (wi > maxDistance) {
                    maxDistance = wi;
                    maxDistance_idx1 = i;
                    maxDistance_idx2 = j;
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
            mstLineSegment.push_back(ridge::PointPair(p, q));
        }
        masb::intList temp;
        temp.resize(cur_size - 1, cur_sheet);
        mstLineSegment_id.insert(mstLineSegment_id.end(), temp.begin(), temp.end());
        
        
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

        std::cout << "start symplify MST as one polyline by maxDistance" << std::endl;
        Vertex maxDistanceStart = vertex(maxDistance_idx1, g_short);
        Vertex maxDistanceEnd = vertex(maxDistance_idx2, g_short);
        std::cout << "longest path -- " << maxDistance
            << " starts at point -- " << indexmap[maxDistanceStart]<<" ("<< maxDistance_idx1<<")"
            << " ends at point -- " << indexmap[maxDistanceEnd] << " (" << maxDistance_idx2 << ")\n";

        dijkstra_shortest_paths(g_short, maxDistanceStart, predecessor_map(&p[0]).distance_map(&d[0]));

        std::vector< graph_traits< Graph >::vertex_descriptor > path;
        graph_traits< Graph >::vertex_descriptor current = maxDistanceEnd;
        
        while (current != maxDistanceStart) {
            path.push_back(current);
            current = p[current];
        }
        path.push_back(maxDistanceStart);
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
        polylines_maxDistance.push_back(a_line);
        polyline_id.push_back(cur_sheet);

        std::cout << "start symplify MST as one polyline by maxAccDist and maxPtNum" << std::endl;
        float maxAccDist = 0;
        int maxPtNum;
        Vertex maxAccDistStart, maxAccDistEnd;
        Vertex maxPtNumStart, maxPtNumEnd;
        for (std::size_t j = 0; j < cur_size; ++j) {
            Vertex s = vertex(j, g_short);
            //dijkstra_shortest_paths(g_short, s, &p[0], &d[0], weightmap_s, indexmap,
            //    std::less<int>(), closed_plus<int>(),
            //    (std::numeric_limits<int>::max)(), 0,
            //    default_dijkstra_visitor());
            dijkstra_shortest_paths(g_short, s, predecessor_map(&p[0]).distance_map(&d[0]));
            //std::cout << "distances and parents:" << std::endl;
            //graph_traits < Graph >::vertex_iterator vi, vend;
            //for (tie(vi, vend) = vertices(g_short); vi != vend; ++vi) {
            //    std::cout << "distance(" << indexmap[*vi] << ") = " << d[*vi] << ", ";
            //    std::cout << "parent(" << indexmap[*vi] << ") = " << indexmap[p[*vi]] << std::
            //        endl;
            //}
            //std::cout << std::endl;

            graph_traits < Graph >::vertex_iterator vii, vendd;
            for (tie(vii, vendd) = vertices(g_short); vii != vendd; ++vii) {
                if (maxAccDist < d[*vii]) {
                    maxAccDist = d[*vii];
                    maxAccDistEnd = *vii;
                    maxAccDistStart = s;
                }
            }
        }
        //because d,p are modified, need to calculate shortest path again
        dijkstra_shortest_paths(g_short, maxAccDistStart, predecessor_map(&p[0]).distance_map(&d[0]));

        std::vector< graph_traits< Graph >::vertex_descriptor > path2;
        graph_traits< Graph >::vertex_descriptor current2 = maxAccDistEnd;

        while (current2 != maxAccDistStart) {
            path2.push_back(current2);
            current2 = p[current2];
        }
        path2.push_back(maxAccDistStart);
        masb::PointList a_line2;
        std::vector< graph_traits< Graph >::vertex_descriptor >::iterator it2;
        for (it2 = path2.begin(); it2 != path2.end(); ++it2) {
            //std::cout << indexmap[*it] << " ";
            auto pt_idx = indexmap[*it2];
            masb::Point pt = cur_candidate[pt_idx];
            //std::cout << pt[0] << " " << pt[1] << " " << pt[2] << std::endl;
            a_line2.push_back(pt);
        }
        //std::cout << std::endl;
        polylines_maxAccDist.push_back(a_line2);

    }
}
/*
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
*/