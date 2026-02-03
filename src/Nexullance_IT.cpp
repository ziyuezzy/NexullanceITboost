#include "Nexullance_IT.hpp"
#include "graph_utility.hpp"
#include <boost/graph/adjacency_list.hpp>
// #include <boost/property_map/property_map.hpp>
#include <iostream>
#include <algorithm>
#include <numeric>
#include <iterator>
#include <random>
// #include <pybind11/pybind11.h>
// #include <pybind11/stl.h>

Nexullance_IT::Nexullance_IT(Graph& _input_graph, const float** _M_EPs, 
            const int _EPR, const float _Cap_core, const float _Cap_access, const bool _verbose): 
            G(_input_graph), Cap_core(_Cap_core), Cap_access(_Cap_access), verbose(_verbose), EPR(_EPR) {

    num_edges = boost::num_edges(G);
    num_vertices = boost::num_vertices(G);

    std::pair<float, float> temp_result=procress_M_EPs(_M_EPs, num_vertices, EPR, &M_R);
    // Calculate the max_access_link_load for later usage
    float maxEPflow = temp_result.first;
    total_flow = temp_result.second;
    // Calculate the max_access_load
    max_access_load = maxEPflow/Cap_access;
    //======================

    next_path_id = 0;
    routing_tables = new std::unordered_map<path_id, float>*[num_vertices];
    for (int i = 0; i < num_vertices; i++) {
        routing_tables[i] = new std::unordered_map<path_id, float>[num_vertices];
    }

    // Initialize link_load and link_path_ids for all edges
    EdgeIterator ei, ei_end;
    for (tie(ei, ei_end) = boost::edges(G); ei != ei_end; ei++) {
        link_load[*ei] = 0.0f;
        link_path_ids[*ei] = std::vector<path_id>();
    }

    // initialize data structures
    all_paths_all_s_d = new std::vector<std::vector<Vertex>>*[num_vertices];
    for (int i = 0; i < num_vertices; i++) {
        all_paths_all_s_d[i] = new std::vector<std::vector<Vertex>>[num_vertices];
    }
    weightmap = get(edge_weight, G);
}

Nexullance_IT::~Nexullance_IT() {
    for (int i = 0; i < num_vertices; i++) {
        delete[] routing_tables[i];
        delete[] all_paths_all_s_d[i];
        delete[] M_R[i];
    }
    delete[] M_R;
    delete[] routing_tables;
    delete[] all_paths_all_s_d;
}

void Nexullance_IT::step_1(float _alpha, float _beta) {
    // first calculate the paths
    
    compute_all_shortest_paths_all_s_d(G, all_paths_all_s_d, weightmap);
    // if(verbose)
    //     std::cout<<"step 1: computed all shortest paths all s d"<<std::endl;

    // then set all loads to 0.0, empty path ids, empty path_id_to_path
    path_id_to_path.clear();
    EdgeIterator ei, ei_end;
    for (tie(ei, ei_end) = boost::edges(G); ei!= ei_end; ei++) {
        link_load[*ei] = 0.0f;
        link_path_ids[*ei].clear();
    }
    next_path_id = 0;

    // if(verbose)
    //     std::cout<<"step 1: cleared link load and path ids"<<std::endl;
    
    // clear and update routing table, update link load, update path ids
    for (int s = 0; s < num_vertices; s++) {
        for (int d = 0; d < num_vertices; d++) {
            routing_tables[s][d].clear();
            if(s==d){
                continue;
            }
            std::vector<std::vector<Vertex>> paths = all_paths_all_s_d[s][d];
            // if(verbose)
            //     std::cout<<"step 1: cleared routing table for s = "<<s<<" d = "<<d<<std::endl;

            float ECMP_weight = 1.0f / paths.size();
            for (const std::vector<Vertex>& path : paths) {
                path_id current_path_id = next_path_id++;
                path_id_to_path[current_path_id] = path;
                routing_tables[s][d][current_path_id] = ECMP_weight;
                for (size_t l = 0; l < path.size() - 1; l++) {
                    Vertex u = path[l];
                    Vertex v = path[l+1];
                    Edge e = boost::edge(u, v, G).first; // the boost::edge(u, v, G) returns a pair<Edge edge(u,v), bool found>
                    link_path_ids[e].push_back(current_path_id);
                    link_load[e] += ECMP_weight*M_R[s][d]/Cap_core; // TODO: to further optimize, divide Cap_core on the M_R at the beginning?
                    boost::put(boost::edge_weight, G, e, _alpha + pow(link_load[e],_beta)); // TODO: alternatively, update the weights outside this loop, may lead to speedup
                }
            }
        }
    }

    // if(verbose)
    //     std::cout<<"step 1: updated routing table, update link load, path ids, and weights"<<std::endl;


    // Find the max value of link_load
    float max_load = compute_max_link_load();
    result_max_loads_step_1.push_back(max_load);
    if(verbose)
        std::cout<<"result from step1: " << result_max_loads_step_1.back() << std::endl;
    return;
}

bool Nexullance_IT::step_2(float _alpha, float _beta, float step, float threshold, int min_attempts, int max_attempts) {
    
    auto rng = std::default_random_engine{};
    int attempts = 0;
    std::list<float> max_loads;

    while (attempts < max_attempts) {
        // Find the max value of link_load
        float max_load = compute_max_link_load();
        max_loads.push_back(max_load);
        result_max_load_step_2=max_load;

        if(max_load <= max_access_load){
            if (verbose){
                std::cout<<"max core link load reaches max access link load, terminating "<<std::endl;
            }
            result_max_load_step_2=max_load;
            return false;
        }

        if (verbose) {
            std::cout << "step2, it=" << attempts << ", max_load = " << max_load 
                      << ", estimated phi = " << total_flow / std::max(max_access_load, result_max_load_step_2) / (num_vertices * EPR) 
                      << std::endl;
        }

        // Check if progress is too slow
        if (attempts > min_attempts) {
            float recent_avg = std::accumulate(std::prev(max_loads.end(), min_attempts / 2), 
                                               max_loads.end(), 0.0f) / static_cast<float>(min_attempts);
            if ((recent_avg - max_load) < threshold) {
                if (verbose) {
                    std::cout << "step 2: low progress, terminating for step = " << step << std::endl;
                    std::cout << "step 2: found max link load = " << max_load << std::endl;
                }
                num_attempts_step_2 += attempts;
                return true;
            }
        }
        
        bool success_attempt = false;
        EdgeIterator ei, ei_end;
        
        for (tie(ei, ei_end) = boost::edges(G); ei != ei_end; ei++) {
            Edge e = *ei;
            if (link_load[e] < max_load) {
                continue;
            }
            const std::vector<path_id>& path_ids = link_path_ids[e];

            std::multimap<float, path_id, std::greater<float>> sorted_path_ids;
            for (const auto& pid : path_ids) {
                const std::vector<Vertex>& path = path_id_to_path[pid];
                Vertex src = path.front();
                Vertex dst = path.back();
                float contribution = routing_tables[src][dst][pid] * M_R[src][dst] / Cap_core;
                sorted_path_ids.insert(std::make_pair(contribution, pid));
            }

            for (const auto& item : sorted_path_ids) {
                path_id old_path_id = item.second;
                const std::vector<Vertex>& old_path = path_id_to_path[old_path_id];
                Vertex src = old_path.front();
                Vertex dst = old_path.back();

                if (verbose) {
                    std::cout << "step 2: starting with old path: ";
                    for (auto v : old_path) {
                        std::cout << v << " ";
                    }
                    std::cout << std::endl;
                }
                        
                std::vector<std::vector<Vertex>> all_paths;
                compute_all_shortest_paths_single_s_d(G, src, dst, all_paths, weightmap);

                if (verbose) {
                    std::cout << "step 2: found " << all_paths.size() << " new paths for s = " << src << " d = " << dst << std::endl;
                }

                std::shuffle(std::begin(all_paths), std::end(all_paths), rng);

                for (const std::vector<Vertex>& new_path : all_paths) {
                    if (new_path != old_path) {
                        success_attempt = true;
                        attempts++;

                        if (verbose) {
                            std::cout << "step 2: starting with new path: ";
                            for (auto v : new_path) {
                                std::cout << v << " ";
                            }
                            std::cout << std::endl;
                        }
                        
                        // Update the paths
                        float delta_weight = 0.0f;
                        std::unordered_map<path_id, float>& current_routing_table = routing_tables[src][dst];

                        auto iter = current_routing_table.find(old_path_id);
                        assert(iter != current_routing_table.end());
                        float old_path_weight = iter->second;

                        bool new_path_found = false;
                        path_id new_path_id = 0;
                        for (auto it = current_routing_table.begin(); it != current_routing_table.end(); ++it) {
                            if (path_id_to_path[it->first] == new_path) {
                                new_path_found = true;
                                new_path_id = it->first;
                                float prev_path_weight = it->second;
                                delta_weight = std::min(step, std::min(old_path_weight, 1.0f - it->second));
                                current_routing_table.erase(it);
                                current_routing_table.insert(std::make_pair(new_path_id, prev_path_weight + delta_weight));
                                break;
                            }
                        }
                        if (!new_path_found) {
                            new_path_id = next_path_id++;
                            path_id_to_path[new_path_id] = new_path;
                            delta_weight = std::min(step, old_path_weight);
                            current_routing_table.insert(std::make_pair(new_path_id, delta_weight));
                        }

                        // Update link loads on the new path
                        for (size_t l = 0; l < new_path.size() - 1; l++) {
                            Vertex u = new_path[l];
                            Vertex v = new_path[l + 1];
                            Edge e = boost::edge(u, v, G).first;
                            link_load[e] += delta_weight * M_R[src][dst] / Cap_core;
                            boost::put(boost::edge_weight, G, e, _alpha + pow(link_load[e], _beta));
                            if (!new_path_found) {
                                link_path_ids[e].push_back(new_path_id);
                            }
                        }
                        
                        current_routing_table.erase(old_path_id);
                        current_routing_table.insert(std::make_pair(old_path_id, old_path_weight - delta_weight));

                        // Update link loads on the old path
                        constexpr float EPSILON = 1e-6f;
                        for (size_t l = 0; l < old_path.size() - 1; l++) {
                            Vertex u = old_path[l];
                            Vertex v = old_path[l + 1];
                            Edge e = boost::edge(u, v, G).first;
                            link_load[e] -= delta_weight * M_R[src][dst] / Cap_core;
                            boost::put(boost::edge_weight, G, e, _alpha + pow(link_load[e], _beta));
                        
                            if (current_routing_table[old_path_id] < EPSILON) {
                                auto iter = std::find(link_path_ids[e].begin(), link_path_ids[e].end(), old_path_id);
                                if (iter != link_path_ids[e].end()) {
                                    link_path_ids[e].erase(iter);
                                }
                            }
                        }
                        
                        if (current_routing_table[old_path_id] < EPSILON) {
                            current_routing_table.erase(old_path_id);
                            path_id_to_path.erase(old_path_id);
                        }
                        
                        break;
                    }
                }

                if (success_attempt) {
                    break;
                }
            }

            if (success_attempt) {
                break;
            }
        }

        if (!success_attempt) {
            if (verbose) {
                std::cout << "step 2: found max link load = " << max_load << std::endl;
                std::cout << "step 2: no progress, terminating after " << attempts << " attempts" << std::endl;
            }
            result_max_load_step_2 = max_load;
            num_attempts_step_2 += attempts;
            return false;
        }
    }

    std::cout << "step 2: max number of attempts reached with step = " << step 
              << ", threshold = " << threshold << ", min_attempts = " << min_attempts 
              << ", max_attempts = " << max_attempts << std::endl;
    num_attempts_step_2 += attempts;
    return true;
}

void Nexullance_IT::optimize(int num_step_1, float alpha_step_1, float beta_step_1, int max_num_step_2, 
        float alpha_step_2, float beta_step_2, int step_2_min_attempts, int step_2_stepping_threshold, int step_2_max_attempts) {
    assert(num_step_1 >= 1 && "num_step_1 should be a positive integer >= 1.");
    
    // Run step 1 multiple times
    for (int i = 0; i < num_step_1; i++) {
        step_1(alpha_step_1, beta_step_1);
    }
    
    // Run step 2 with decreasing step sizes
    float step = 0.5f;
    for (int i = 0; i < max_num_step_2; i++) {
        if (step_2(alpha_step_2, beta_step_2, step, step_2_stepping_threshold, step_2_min_attempts, step_2_max_attempts)) {
            step *= 0.5f;
        } else {
            break;
        }
    }
}


result_routing_table Nexullance_IT::get_routing_table() {
    result_routing_table result;
    
    // Iterate over routing_tables and convert to result_routing_table
    for (size_t s = 0; s < num_vertices; s++) {
        for (size_t d = 0; d < num_vertices; d++) {
            if (s == d) {
                continue;
            }

            std::vector<std::pair<std::vector<Vertex>, float>> paths;
            for (const auto& item : routing_tables[s][d]) {
                const std::vector<Vertex>& path = path_id_to_path.at(item.first);
                paths.emplace_back(path, item.second);
            }
            result.insert(std::make_pair(std::make_pair(s,d), paths));
        }
    }
    return result;
}

float Nexullance_IT::get_max_core_load() const {
    return result_max_load_step_2;
}

float Nexullance_IT::get_phi() const {
    return total_flow / std::max(max_access_load, result_max_load_step_2) / (num_vertices * EPR);
}

float Nexullance_IT::compute_max_link_load() const {
    float max_load = 0.0f;
    EdgeIterator ei, ei_end;
    for (tie(ei, ei_end) = boost::edges(G); ei != ei_end; ei++) {
        max_load = std::max(max_load, link_load.at(*ei));
    }
    return max_load;
}

void Nexullance_IT::update_edge_weight(Edge e, float alpha, float beta) {
    boost::put(boost::edge_weight, G, e, alpha + pow(link_load[e], beta));
}
