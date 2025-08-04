// actually this includes SD_Nexullance_IT and diff_Nexullance_IT,
// the difference is only in a flag, to control which starting routing table to use.

#include "Nexullance_IT_fast.hpp"
#include "graph_utility.hpp"
#include <chrono>
#include <algorithm>
#include <functional>
#include <numeric>
#include <iostream>
#include <random>

// Constructor implementation
Nexullance_IT_fast::Nexullance_IT_fast(Graph& _input_graph, const float _Cap_core,
                                     const float _Cap_access, const bool _verbose) 
    : G(_input_graph), Cap_core(_Cap_core), Cap_access(_Cap_access), verbose(_verbose) {
    num_edges = boost::num_edges(G);
    num_vertices = boost::num_vertices(G);

    // Initialize link_load for all edges
    EdgeIterator ei, ei_end;
    for (boost::tie(ei, ei_end) = boost::edges(G); ei != ei_end; ++ei) {
        Edge e = *ei;
        link_load.emplace(e, 0.0f);
    }

    rng = std::default_random_engine {};
}

// Destructor implementation
Nexullance_IT_fast::~Nexullance_IT_fast() {
    // Clean up dynamically allocated memory for path_info objects
    for (auto& src_vec : initial_RT) {
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
void Nexullance_IT_fast::initialization(size_t max_path_length, size_t initial_weight_max_path_length, bool disjoint_paths) {
    if (verbose) {
        std::cout << "Initializing with max path length: " << max_path_length 
                  << ", initial weight max path length: " << initial_weight_max_path_length << std::endl;
    }

    assert(max_path_length > 2 && "Max path length must be greater than 2.");
    
    initial_RT.resize(num_vertices);
    for (size_t i = 0; i < num_vertices; ++i)
        initial_RT[i].resize(num_vertices);

    std::vector<std::vector<std::vector<std::vector<Vertex>>>> all_paths_matrix;

    // Find all paths for all source-destination pairs up to max_path_length
    std::vector<std::vector<std::vector<std::vector<Vertex>>>> all_paths_matrix_unfiltered;
    compute_all_paths_with_max_length_all_s_d(G, max_path_length, all_paths_matrix_unfiltered);


    if (disjoint_paths){
        std::vector<std::vector<std::vector<std::vector<Vertex>>>> all_paths_matrix_filtered;
        // Filter paths to be edge-disjoint
        std::vector<std::vector<std::vector<std::vector<Vertex>>>> filtered_paths;
        filter_edge_disjoint_paths(all_paths_matrix_unfiltered, all_paths_matrix_filtered, G, true);
        all_paths_matrix = std::move(all_paths_matrix_filtered);
    }else{
        all_paths_matrix = std::move(all_paths_matrix_unfiltered);
    }
    
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
                    initial_RT[src][dst].push_back(new_path_info);
                    
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
                    initial_RT[src][dst].push_back(new_path_info);
                    
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

    previous_result_RT = deep_copy_RT(initial_RT); // Set previous result as a copy of the initial routing table
    current_RT = deep_copy_RT(initial_RT); // Set current result as a copy of the initial routing table
    
    // Initialize Map_passby for all edges
    EdgeIterator ei, ei_end;
    for (boost::tie(ei, ei_end) = boost::edges(G); ei != ei_end; ++ei) {
        Edge e = *ei;
        Map_passby.emplace(e, descend_loadcontri_pathinfos());
        
        std::vector<std::vector<std::vector<path_info*>>> not_passby_3d;
        not_passby_3d.resize(num_vertices);
        for (size_t src = 0; src < num_vertices; ++src) {
            not_passby_3d[src].resize(num_vertices);
        }
        Map_NOT_passby.emplace(e, std::move(not_passby_3d));
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
            for (auto& path_ptr : current_RT[src][dst]) {
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
                    auto passby_it = Map_passby.find(found_edge);
                    #ifdef DEBUG
                    assert(passby_it != Map_passby.end() && "Edge not found in Map_passby.");
                    #endif
                    passby_it->second.insert({path_ptr, 0.0f});

                    path_edges.insert(found_edge);
                    
                    ++it;
                    ++next_it;
                }
                
                // For edges not in this path, add to Map_NOT_passby
                for (const Edge& e : all_edges) {
                    if (path_edges.find(e) == path_edges.end()) {
                        auto not_passby_it = Map_NOT_passby.find(e);
                        #ifdef DEBUG
                        assert(not_passby_it != Map_NOT_passby.end() && "Edge not found in Map_NOT_passby.");
                        #endif
                        not_passby_it->second[src][dst].push_back(path_ptr);
                    }
                }
            }
        }
    }
    
    if (verbose) {
        std::cout << "Initialization complete." << std::endl;
    }
}

// Optimize for M_EPs method implementation
IT_outputs Nexullance_IT_fast::optimize_for_M_EPs(float** M_EPs, size_t EPR, size_t max_num_fix_step, bool from_initial_RT) {
    // Start timing
    auto start_time = std::chrono::high_resolution_clock::now();

    // Process the M_EPs matrix to get M_R
    float** M_R = nullptr;
    auto temp_result = process_M_EPs(const_cast<const float**>(M_EPs), num_vertices, EPR, &M_R);
    float max_EP_flow = temp_result.first;
    float total_flow = temp_result.second;
    float max_access_load = max_EP_flow / Cap_access;

    // set starting routing table, either from scratch (original nexu) from the previous result (diff nexu)
    if (from_initial_RT)
        set_RT_as_RT(current_RT, initial_RT); // Use initial routing table
    else 
        set_RT_as_RT(current_RT, previous_result_RT); // Use previous result routing table


    // first set all link load to zero
    for (auto& load_pair : link_load) {
        load_pair.second = 0.0f; // Reset all link loads to zero
    }

    // calculate based on the new input matrix
    for (int s = 0; s < num_vertices; s++) {
        for (int d = 0; d < num_vertices; d++) {
            for (path_info* path_info_it: current_RT[s][d]) {
                path current_path = path_info_it->first;
                float path_weight = path_info_it->second;
                float path_load_contri = path_weight * M_R[s][d] / Cap_access; // the contribution of the flow on this path to the load on the links
                
                auto it_u = current_path.begin();
                for (int l = 0; l < current_path.size() - 1; l++) {
                    auto it_v = it_u;
                    std::advance(it_v, 1);
                    
                    Vertex u = *it_u;
                    Vertex v = *it_v;
                    Edge e = boost::edge(u, v, G).first; // the boost::edge(u, v, G) returns a pair<Edge edge(u,v), bool found>
                    
                    // modify path load contribution, if higher than 0.0f
                    if (path_load_contri >= 0.0f) {
                        auto Map_passby_it = Map_passby.find(e);
                        #ifdef DEBUG
                        assert(Map_passby_it != Map_passby.end() && "Edge not found in Map_passby.");
                        #endif
                        auto& key_index = Map_passby_it->second.get<by_key>();
                        auto key_it = key_index.find(path_info_it);
                        #ifdef DEBUG
                        assert(key_it != key_index.end() && "Path not found in Map_passby for this edge.");
                        #endif
                        // Update the load contribution for this path in Map_passby
                        key_index.modify(key_it, [path_load_contri](Entry& entry) {
                            entry.value = path_load_contri;
                        });

                    }
                    
                    // calculate link load
                    auto load_it = link_load.find(e);
                    #ifdef DEBUG
                    assert(load_it != link_load.end() && "Edge not found in link_load.");
                    #endif
                    load_it->second += path_load_contri;

                    it_u = it_v;
                }
            }
        }
    }

    // do the steps
    bool to_continue = true;
    size_t tot_attempts = 0;
    size_t num_attempts;
    float max_core_load;
    float step = 0.5;
    size_t num_fix_step = 0;
    while ((max_num_fix_step > num_fix_step) && to_continue)
    {
        std::tie(to_continue, num_attempts, max_core_load) = optimize_for_M_R_fixed_step(M_R, 
            step, max_access_load, total_flow);
        step *= 0.5;
        tot_attempts += num_attempts;
        num_fix_step++;
    }

    // end timer
    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end_time - start_time;
    if(verbose)     std::cout<<"elapsed time = " << elapsed.count() << " s" << std::endl;
    float phi = total_flow/std::max(max_access_load, max_core_load)/(num_vertices*EPR);
    
    // Cleanup M_R
    for (size_t i = 0; i < num_vertices; ++i)   delete[] M_R[i];
    delete[] M_R;
    set_RT_as_RT(previous_result_RT, current_RT); // Set previous result to current routing table

    IT_outputs result = IT_outputs(elapsed.count(), max_core_load, phi, get_routing_table(previous_result_RT), tot_attempts);
    return result;
}

// return <bool continue or not, num_attempts, max_core_load>
std::tuple<bool, size_t, float> Nexullance_IT_fast::optimize_for_M_R_fixed_step(float** M_R, float step, float max_access_load, float total_flow){
    
    if (verbose)
        std::cout << "Starting optimization with step = " << step << std::endl;

    // start the main loop
    size_t attempts = 0;
    std::list<float> max_loads_hist;
    float result_max_load = 0.0;
    // path_info* last_decreased_path = nullptr;
    // path_info* last_increased_path = nullptr;
    
    while (attempts < max_attempts) {
        // Find max load link
        float max_load = 0.0;
        std::vector<Edge> max_load_edges;
        
        EdgeIterator ei, ei_end;
        for (boost::tie(ei, ei_end) = boost::edges(G); ei != ei_end; ++ei) {
            Edge e = *ei;
            auto load_it = link_load.find(e);
            #ifdef DEBUG
            assert(load_it != link_load.end() && "Edge not found in link_load.");
            #endif
            if (load_it->second > max_load) {
                max_load = load_it->second;
                max_load_edges.clear();
                max_load_edges.push_back(e);
            } else if (load_it->second == max_load) {
                max_load_edges.push_back(e);
            }
        }
        
        max_loads_hist.push_back(max_load);
        result_max_load = max_load;
        
        if (max_load <= max_access_load) {
            if (verbose) {
                std::cout << "Max core link load reaches max access link load, terminating" << std::endl;
            }
            return std::make_tuple(false, attempts, result_max_load);
        }
        
        if (verbose) {
            std::cout << "Step = " << step << ", iteration = " << attempts 
                      << ", max_load = " << result_max_load 
                    //   << ", estimated phi = " << total_flow / std::max(max_access_load, result_max_load) / (num_vertices * EPR) 
                      << std::endl;
        }
        
        // Check if we've made sufficient progress
        if(attempts > min_attempts){
            float recent_average_load = std::accumulate(std::prev(max_loads_hist.end(), min_attempts/2), max_loads_hist.end(), 0.0f)/((float)(min_attempts/2));
            if  ((recent_average_load - max_load)<progress_threshold){
                if (verbose){
                    std::cout<<"nexu fast: low progress, terminating for step = "<< step <<std::endl;
                    std::cout<<"nexu fast: found max link load" << result_max_load <<std::endl;
                }
                return std::make_tuple(true, attempts, result_max_load);
            }
        }
        
        bool success_attempt = false;
        
        // Shuffle max_load_edges to try them in random order
        std::shuffle(max_load_edges.begin(), max_load_edges.end(), rng);
        
        for (const Edge& max_edge : max_load_edges) {
            // Get all paths passing by this edge and sort them by contribution to load
            auto passby_it = Map_passby.find(max_edge);
            #ifdef DEBUG
            assert(passby_it != Map_passby.end() && "Edge not found in Map_passby.");
            assert(!passby_it->second.empty() && "No paths passing through this edge, which is not possible.");
            #endif
            
            // TODO: abstract the path-sorting logic into a separate function here.
                // 1. the current implementation is sorting paths by their contribution to load
                // 2. another implementation can just shuffle the paths randomly
            
            // std::vector<path_info*>& paths_through_edge = passby_it->second;
            auto& sorted_index = passby_it->second.get<by_value>();

            // // Create a multimap to sort paths by their contribution to load
            // std::multimap<float, path_info*, std::greater<float>> sorted_paths;
            // for (path_info* path_ptr : paths_through_edge) {
            //     // Find source and destination
            //     Vertex src = path_ptr->first.front();
            //     Vertex dst = path_ptr->first.back();
            //     // Calculate contribution to load
            //     float contribution = path_ptr->second * M_R[src][dst] / Cap_core;
            //     sorted_paths.insert(std::make_pair(contribution, path_ptr));
            // }

            if (verbose) {
                std::cout << "Max load edge: " << max_edge 
                          << ", number of paths through this edge: " << sorted_index.size() << std::endl;
            }
            // Try to reroute each path, starting with those contributing most to load
            for (auto& pathinfo_and_contri : sorted_index) {
                path_info* old_path_ptr = pathinfo_and_contri.key;
                Vertex src = old_path_ptr->first.front();
                Vertex dst = old_path_ptr->first.back();

                // TODO: can use different strategies to select an alternative path.
                auto select_result = select_alternative_path_random(old_path_ptr, max_edge, max_load);
                path_info* new_path_ptr = select_result.first;
                float new_path_max_load = select_result.second;
                if (new_path_ptr == nullptr)
                    continue;  // no alternative path will improve the situation

                // Calculate how much weight to transfer
                float delta_weight = std::min(std::min(step, old_path_ptr->second), 1-new_path_ptr->second);
                delta_weight = std::min(delta_weight, Cap_core * (max_load - new_path_max_load) / M_R[src][dst]);
                
                if (verbose)
                {
                    std::cout << "Attempting to reroute path from " << src << " to " << dst 
                                << " with old path: " << path_to_string(old_path_ptr->first)
                                << ", new path: " << path_to_string(new_path_ptr->first)
                              << ", old path max load: " << max_load 
                              << ", new path max load: " << new_path_max_load 
                              << ", old path weight: " << old_path_ptr->second 
                              << ", new path weight: " << new_path_ptr->second << std::endl;
                }
                
                // Update path weights
                old_path_ptr->second -= delta_weight;
                new_path_ptr->second += delta_weight;

                // calculate the new load contributions
                float old_path_load_contri = old_path_ptr->second * M_R[src][dst] / Cap_core;
                float new_path_load_contri = new_path_ptr->second * M_R[src][dst] / Cap_core;

                // Update link loads and load contributions in Map_passby
                // Go through the old path
                auto it_old_u = old_path_ptr->first.begin();
                for (int l = 0; l < old_path_ptr->first.size() - 1; l++) {
                    auto it_old_v = it_old_u;
                    std::advance(it_old_v, 1);
                    
                    Vertex u = *it_old_u;
                    Vertex v = *it_old_v;
                    Edge e = boost::edge(u, v, G).first;
                    
                    // Update link load
                    auto load_it = link_load.find(e);
                    #ifdef DEBUG
                    assert(load_it != link_load.end() && "Edge not found in link_load.");
                    #endif
                    load_it->second -= delta_weight * M_R[src][dst] / Cap_core;
                    // handle rounding errors
                    #ifdef DEBUG
                    assert(load_it->second > -1E-6f && "Link load should not be negative");
                    #endif
                    if (load_it->second < 1E-6f) load_it->second = 0.0f;

                    // Update the load contribution in Map_passby
                    auto Map_passby_it = Map_passby.find(e);
                    #ifdef DEBUG
                    assert(Map_passby_it != Map_passby.end() && "Edge not found in Map_passby.");
                    #endif
                    auto& key_index = Map_passby_it->second.get<by_key>();
                    auto key_it = key_index.find(old_path_ptr);
                    #ifdef DEBUG
                    assert(key_it != key_index.end() && "Path not found in Map_passby for this edge.");
                    #endif
                    // Update the load contribution for this path in Map_passby
                    key_index.modify(key_it, [old_path_load_contri](Entry& entry) {
                        entry.value = old_path_load_contri;
                    });

                    it_old_u = it_old_v;
                }

                // Go through the new path
                auto it_new_u = new_path_ptr->first.begin();
                for (int l = 0; l < new_path_ptr->first.size() - 1; l++) {
                    auto it_new_v = it_new_u;
                    std::advance(it_new_v, 1);
                    
                    Vertex u = *it_new_u;
                    Vertex v = *it_new_v;
                    Edge e = boost::edge(u, v, G).first;

                    // update link load                    
                    auto load_it = link_load.find(e);
                    #ifdef DEBUG
                    assert(load_it != link_load.end() && "Edge not found in link_load.");
                    #endif
                    load_it->second += delta_weight * M_R[src][dst] / Cap_core;

                    // update the load contribution in Map_passby
                    auto Map_passby_it = Map_passby.find(e);
                    #ifdef DEBUG
                    assert(Map_passby_it != Map_passby.end() && "Edge not found in Map_passby.");
                    #endif
                    auto& key_index = Map_passby_it->second.get<by_key>();
                    auto key_it = key_index.find(new_path_ptr);
                    #ifdef DEBUG
                    assert(key_it != key_index.end() && "Path not found in Map_passby for this edge.");
                    #endif
                    // Update the load contribution for this path in Map_passby
                    key_index.modify(key_it, [new_path_load_contri](Entry& entry) {
                        entry.value = new_path_load_contri;
                    });
                    
                    it_new_u = it_new_v;
                }
                success_attempt = true;
                attempts++;
                
                if (success_attempt) {
                    break;  // Found and applied a change, move to next iteration
                }
            }
            
            if (success_attempt) {
                break;  // Found and applied a change, move to next iteration
            }
        }
        
        if (!success_attempt) {
            if (verbose) {
                std::cout << "No progress possible, terminating after " << attempts 
                          << " attempts with step = " << step << std::endl;
                std::cout << "Found max link load " << result_max_load << std::endl;
            }
            return std::make_tuple(false, attempts, result_max_load);
        }
    }
    
    if (verbose) {
        std::cout << "Max number of attempts reached with step = " << step << std::endl;
    }
    
    return std::make_tuple(true, attempts, result_max_load);
}


std::pair<path_info*, float> Nexullance_IT_fast::select_alternative_path_random(path_info* old_path_info, Edge max_load_edge, float old_path_max_load){

    // Get paths that don't pass through this edge for this src-dst pair
    auto not_passby_it = Map_NOT_passby.find(max_load_edge);
    #ifdef DEBUG
    assert(not_passby_it != Map_NOT_passby.end() && "Edge not found in Map_NOT_passby.");
    #endif
    
    Vertex src = old_path_info->first.front();
    Vertex dst = old_path_info->first.back();
    
    std::vector<path_info*>& alternative_paths = not_passby_it->second[src][dst];
    
    if (alternative_paths.empty()) {
        return {nullptr, 0.0f};  // No alternative paths, continue the loop
    }

    // Shuffle alternative paths to try them in random order
    std::shuffle(alternative_paths.begin(), alternative_paths.end(), rng);
    
    for (path_info* new_path_ptr : alternative_paths) {
        // // Avoid oscillation by checking if we're reversing a recent change
        // if ((old_path_ptr == last_increased_path) && (new_path_ptr == last_decreased_path)) {
        //     continue;
        // }
        
        // Find the max load along the new path
        float new_path_max_load = 0.0f;
        auto it_u = new_path_ptr->first.begin();
        for (int l = 0; l < new_path_ptr->first.size() - 1; l++) {
            auto it_v = it_u;
            std::advance(it_v, 1);
            
            Vertex u = *it_u;
            Vertex v = *it_v;
            Edge e = boost::edge(u, v, G).first;
            
            auto load_it = link_load.find(e);
            #ifdef DEBUG
            assert(load_it != link_load.end() && "Edge not found in link_load.");
            #endif
            if (load_it->second > new_path_max_load) {
                new_path_max_load = load_it->second;
            }
            
            it_u = it_v;
        }
        
        #ifdef DEBUG
        assert(new_path_max_load <= old_path_max_load && "New path max load should not exceed current max load");
        #endif
        if (new_path_max_load == old_path_max_load) {
            continue;  // New path doesn't help, try another one
        }else{
            return {new_path_ptr, new_path_max_load};  // Found a suitable alternative path
        }
    }
    return {nullptr, 0.0f};  // No suitable alternative path found, return empty pair
}


result_routing_table Nexullance_IT_fast::get_routing_table(
    std::vector< std::vector< std::vector< path_info* > > > routing_table){
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

    ram_usage += sizeof(previous_result_RT);
    ram_usage += sizeof(Map_passby);
    ram_usage += sizeof(Map_NOT_passby);

    // Add size of path_info variables
    for (const auto& src_vec : previous_result_RT) {
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
