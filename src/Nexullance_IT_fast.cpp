// actually this includes SD_Nexullance_IT and diff_Nexullance_IT,
// the difference is only in a flag, to control which starting routing table to use.

#include "Nexullance_IT_fast.hpp"
#include "graph_utility.hpp"
#include <chrono>
#include <algorithm>
#include <functional>
#include <numeric>
#include <iostream>

// Constructor implementation
Nexullance_IT_fast::Nexullance_IT_fast(Graph& _input_graph, const float _Cap_core,
                                     const float _Cap_access, const bool _verbose) 
    : G(_input_graph), Cap_core(_Cap_core), Cap_access(_Cap_access), verbose(_verbose) {
    num_edges = boost::num_edges(G);
    num_vertices = boost::num_vertices(G);
    
    // Initialize routing table with the correct size
    routing_table.resize(num_vertices);
    for (size_t i = 0; i < num_vertices; ++i)
        routing_table[i].resize(num_vertices);
}

// Destructor implementation
Nexullance_IT_fast::~Nexullance_IT_fast() {
    // Clean up dynamically allocated memory for path_info objects
    for (auto& src_vec : routing_table) {
        for (auto& dst_vec : src_vec) {
            for (auto& path_ptr : dst_vec) {
                delete path_ptr;
            }
        }
    }
    
    if (verbose)
        std::cout << "Nexullance_IT_fast instance destroyed." << std::endl;
}

// Initialization method implementation
// If the second input argument initial_weight_max_path_length is 0, only shortest paths get non-zero initial weights
void Nexullance_IT_fast::initialization(size_t max_path_length, size_t initial_weight_max_path_length) {
    if (verbose) {
        std::cout << "Initializing with max path length: " << max_path_length 
                  << ", initial weight max path length: " << initial_weight_max_path_length << std::endl;
    }

    assert(max_path_length > 2 && "Max path length must be greater than 2.");
    
    // Find all paths for all source-destination pairs up to max_path_length
    std::vector<std::vector<std::vector<std::vector<Vertex>>>> all_paths_matrix;
    compute_all_paths_with_max_length_all_s_d(G, max_path_length, all_paths_matrix);
    
    if (verbose) {
        std::cout << "Finished computing all paths for all source-destination pairs." << std::endl;
    }
    
    // For each source-destination pair, process the pre-calculated paths
    for (size_t src = 0; src < num_vertices; ++src) {
        for (size_t dst = 0; dst < num_vertices; ++dst) {
            if (src == dst) continue; // Skip self-loops
            
            // Get the pre-calculated paths for this source-destination pair
            const auto& all_paths = all_paths_matrix[src][dst];
            
            // Filter paths by length and add to routing table
            std::vector<path_info*> valid_paths;
            
            if (initial_weight_max_path_length == 0 && !all_paths.empty()) {
            // If initial_weight_max_path_length is 0, only shortest paths get non-zero initial weights
                // Find minimum path length
                size_t min_path_length = std::numeric_limits<size_t>::max();
                for (const auto& _path : all_paths) {
                    min_path_length = std::min(min_path_length, _path.size());
                }
                
                // Add only shortest paths to valid_paths
                for (const std::vector<Vertex>& _path : all_paths) {
                    path_info* new_path_info = new path_info(path(_path.begin(), _path.end()), 0.0f);
                    
                    // Add to routing table
                    routing_table[src][dst].push_back(new_path_info);
                    
                    // Add to valid paths only if it's a shortest path
                    if (_path.size() == min_path_length) {
                        valid_paths.push_back(new_path_info);
                    }
                }
            } else {
            // Otherwise, all paths up to initial_weight_max_path_length get non-zero initial weights
                for (const std::vector<Vertex>& _path : all_paths) {
                    path_info* new_path_info = new path_info(path(_path.begin(), _path.end()), 0.0f);
                    
                    // Add to routing table
                    routing_table[src][dst].push_back(new_path_info);
                    
                    // Add to valid paths for initial weight calculation
                    if (_path.size() - 1 <= initial_weight_max_path_length) {
                        valid_paths.push_back(new_path_info);
                    }
                }
            }
            
            // Set initial weights for valid paths
            assert(!valid_paths.empty() && "No valid paths found for source-destination pair.");
            float weight = 1.0f / valid_paths.size();
            for (auto& path_ptr : valid_paths) {
                path_ptr->second = weight;
            }
        }
    }
    
    // Initialize Map_passby for all edges
    EdgeIterator ei, ei_end;
    for (boost::tie(ei, ei_end) = boost::edges(G); ei != ei_end; ++ei) {
        Edge e = *ei;
        Map_passby[e] = std::vector<path_info*>();
        Map_NOT_passby[e].resize(num_vertices);
        for (size_t src = 0; src < num_vertices; ++src) {
            Map_NOT_passby[e][src].resize(num_vertices);
        }
    }
    
    // Process paths once to build both mappings efficiently
    for (size_t src = 0; src < num_vertices; ++src) {
        for (size_t dst = 0; dst < num_vertices; ++dst) {
            if (src == dst) continue;
            
            // Store all edges in the graph for fast checking
            std::set<Edge> all_edges;
            for (boost::tie(ei, ei_end) = boost::edges(G); ei != ei_end; ++ei) {
                all_edges.insert(*ei);
            }
            
            // For each path between src and dst
            for (auto& path_ptr : routing_table[src][dst]) {
                // Get edges in this path
                const path& p = path_ptr->first;
                std::set<Edge> path_edges;
                
                auto it = p.begin();
                auto next_it = std::next(it);
                
                while (next_it != p.end()) {
                    Vertex u = *it;
                    Vertex v = *next_it;
                    
                    Edge found_edge;
                    bool exists;
                    boost::tie(found_edge, exists) = boost::edge(u, v, G);
                    
                    #ifdef DEBUG
                    assert(exists && "Edge does not exist in the graph.");
                    #endif
                    
                    // Add path to Map_passby for this edge
                    Map_passby[found_edge].push_back(path_ptr);
                    path_edges.insert(found_edge);
                    
                    ++it;
                    ++next_it;
                }
                
                // For edges not in this path, add to Map_NOT_passby
                for (const Edge& e : all_edges) {
                    if (path_edges.find(e) == path_edges.end()) {
                        Map_NOT_passby[e][src][dst].push_back(path_ptr);
                    }
                }
            }
        }
    }
    
    if (verbose) {
        std::cout << "Initialization complete." << std::endl;
    }
}

// // Optimize for M_EPs method implementation
// IT_outputs Nexullance_IT_fast::optimize_for_M_EPs(float** M_EPs, float threshold, 
//      size_t max_num_step, int min_attempts, int max_attempts) {
//     // Start timing
//     auto start_time = std::chrono::high_resolution_clock::now();
    
//     // Process the M_EPs matrix to get M_R
//     float** M_R = nullptr;
//     auto result = procress_M_EPs(const_cast<const float**>(M_EPs), num_vertices, EPR, &M_R);
//     float max_EP_flow = result.first;
//     float total_flow = result.second;
    
//     float max_access_load = max_EP_flow / Cap_access;
    
//     if (verbose) {
//         std::cout << "Max EP flow: " << max_EP_flow << ", total flow: " << total_flow 
//                   << ", max access load: " << max_access_load << std::endl;
//     }
    
//     // TODO: fix
//     auto optimization_result = optimize_for_M_R_fixed_step(
//         M_R, threshold, 1.0f / max_num_step, min_attempts, max_attempts, max_access_load, total_flow);
    
//     bool continue_optimization = std::get<0>(optimization_result);
//     size_t num_attempts = std::get<1>(optimization_result);
//     float max_core_link_load = std::get<2>(optimization_result);
    
//     // Calculate phi (throughput)
//     float phi;
//     if (max_access_load > max_core_link_load) {
//         phi = total_flow / max_access_load;
//     } else {
//         phi = total_flow / max_core_link_load;
//     }
    
//     // Stop timing
//     auto end_time = std::chrono::high_resolution_clock::now();
//     std::chrono::duration<double> elapsed = end_time - start_time;
    
//     // Create routing table result. //TODO: check if this is correct
//     result_routing_table rt;
//     for (size_t src = 0; src < num_vertices; ++src) {
//         for (size_t dst = 0; dst < num_vertices; ++dst) {
//             if (src == dst) continue;
            
//             std::pair<Vertex, Vertex> sd_pair = std::make_pair(src, dst);
//             rt[sd_pair] = std::vector<std::pair<std::vector<Vertex>, float>>();
            
//             // Convert path_info to the required format
//             for (auto& path_ptr : routing_table[src][dst]) {
//                 const path& p = path_ptr->first;
//                 float weight = path_ptr->second;
                
//                 std::vector<Vertex> path_vector(p.begin(), p.end());
//                 rt[sd_pair].push_back(std::make_pair(path_vector, weight));
//             }
//         }
//     }
    
//     // Cleanup M_R
//     for (size_t i = 0; i < num_vertices; ++i) {
//         delete[] M_R[i];
//     }
//     delete[] M_R;
    
//     return IT_outputs(elapsed.count(), max_core_link_load, phi, rt, num_attempts);
// }

// // Optimize for M_R fixed step method implementation
// std::tuple<bool, size_t, float> Nexullance_IT_fast::optimize_for_M_R_fixed_step(
//     float** M_R, float threshold, float step, int min_attempts, int max_attempts, 
//     float max_access_load, float total_flow) {
    
//     if (verbose) {
//         std::cout << "Starting optimization with threshold: " << threshold 
//                   << ", step: " << step << std::endl;
//     }

    
//     // return std::make_tuple(continue_optimization, num_attempts, max_core_link_load);
// }


result_routing_table Nexullance_IT_fast::get_routing_table(){
    result_routing_table result = result_routing_table();
    // iterate over "routing_table" and convert it to "result_routing_table"
    for (int s = 0; s < num_vertices; s++) {
        for (int d = 0; d < num_vertices; d++) {
            if (s==d)
                continue;
            
            // Create the source-destination pair key
            std::pair<Vertex, Vertex> sd_pair = std::make_pair(s, d);
            
            // Initialize the vector for this source-destination pair
            result[sd_pair] = std::vector<std::pair<std::vector<Vertex>, float>>();
            
            // Convert each path_info to the required format
            for (auto& path_ptr : routing_table[s][d]) {
                // Extract the path (std::list<Vertex>) and weight from path_info
                const path& p = path_ptr->first;
                float weight = path_ptr->second;
                
                // Convert std::list<Vertex> to std::vector<Vertex>
                std::vector<Vertex> path_vector(p.begin(), p.end());
                
                // Add the path-weight pair to the result
                result[sd_pair].push_back(std::make_pair(path_vector, weight));
            }
        }
    }
    return result;
}

ulong Nexullance_IT_fast::get_RAM_after_init() {
    // Calculate the RAM usage after initialization
    ulong ram_usage = 0; 

    ram_usage += sizeof(routing_table);
    ram_usage += sizeof(Map_passby);
    ram_usage += sizeof(Map_NOT_passby);

    // Add size of path_info variables
    for (const auto& src_vec : routing_table) {
        for (const auto& dst_vec : src_vec) {
            for (const auto& path_ptr : dst_vec) {
                ram_usage += sizeof(*path_ptr) + path_ptr->first.size() * sizeof(Vertex) + sizeof(float);
            }
        }
    }

    return ram_usage;
}
    
    // // Add size of Map_passby and Map_NOT_passby
    // ram_usage += Map_passby.size() * (sizeof(Edge) + sizeof(std::vector<path_info*>));

    // for (const auto& entry : Map_NOT_passby) {
    //     ram_usage += sizeof(Edge);
    //     for (size_t src = 0; src < num_vertices; ++src) {
    //         for (size_t dst = 0; dst < num_vertices; ++dst) {
    //             if (src == dst) continue;
    //             ram_usage += sizeof(std::vector<path_info*>) + entry.second[src][dst].size() * sizeof(path_info*);
    //         }
    //     }
    // }
    