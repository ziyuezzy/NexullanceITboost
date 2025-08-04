#ifndef GRAPH_UTILITY_HPP
#define GRAPH_UTILITY_HPP

#include "definitions.hpp"
#include <Eigen/Dense>
#include <utility>

Graph read_graph_from_arcs(int V, Eigen::MatrixX2i arcs, bool print);

// float compute_network_data_throughput(const float** M_R, float max_core_load);

std::pair<float, float> process_M_EPs(const float** _M_EPs, const int _V, const int _EPR, float*** out_M_R);

// Path-finding algorithm ultilities:
// TODO: a path can be simlified as "Vertex*"
void compute_all_shortest_paths_all_s_d(const Graph &G, std::vector<std::vector<Vertex>>** all_s_d_paths,
                                const boost::property_map<Graph, boost::edge_weight_t>::type &weightmap);

void compute_all_shortest_paths_single_source(const Graph &G, Vertex s, std::vector<std::vector<Vertex>>* all_paths,
                                const boost::property_map<Graph, boost::edge_weight_t>::type &weightmap);

void compute_all_shortest_paths_single_s_d(const Graph &G, Vertex s, Vertex d, std::vector<std::vector<Vertex>> &all_paths,
                                const boost::property_map<Graph, boost::edge_weight_t>::type &weightmap);

// Find all paths with max length for all source-destination pairs in the graph
// // TODO: validate against networkx
void compute_all_paths_with_max_length_all_s_d(const Graph &G, size_t max_length, 
                                             std::vector<std::vector<std::vector<std::vector<Vertex>>>> &all_paths_matrix);

void filter_edge_disjoint_paths(const std::vector<std::vector<std::vector<std::vector<Vertex>>>> &all_paths,
                               std::vector<std::vector<std::vector<std::vector<Vertex>>>> &filtered_paths,
                               const Graph &g,
                               bool verbose = false);

#endif