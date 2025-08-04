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

    std::pair<float, float> temp_result=process_M_EPs(_M_EPs, num_vertices, EPR, &M_R);
    // Calculate the max_access_link_load for later usage
    float maxEPflow = temp_result.first;
    total_flow = temp_result.second;
    // Calculate the max_access_load
    max_access_load = maxEPflow/Cap_access;
    //======================

    next_path_id = 0;
    routing_table = new std::unordered_map<path_id, float>*[num_vertices];
    for (int i = 0; i < num_vertices; i++) {
        routing_table[i] = new std::unordered_map<path_id, float>[num_vertices];
    }

    link_load = new float*[num_vertices];
    for (int i = 0; i < num_vertices; i++) {
        link_load[i] = new float[num_vertices]{0.0f};
    }
    // // assign 0.0 to all link_load initially
    // for (int i = 0; i < num_vertices; i++) { // order of loops??
    //     for (int j = 0; j < num_vertices; j++) {
    //         link_load[i][j] = 0.0;
    //     }
    // }

    link_path_ids = new std::vector<path_id>*[num_vertices];
    for (int i = 0; i < num_vertices; i++) {
        link_path_ids[i] = new std::vector<path_id>[num_vertices];
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
        delete[] routing_table[i];
        delete[] link_load[i];
        delete[] link_path_ids[i];
        delete[] all_paths_all_s_d[i];
        delete[] M_R[i];
    }
    delete[] M_R;
    delete[] routing_table;
    delete[] link_load;
    delete[] link_path_ids;
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
        link_load[(*ei).m_source][(*ei).m_target] = 0;
        link_path_ids[(*ei).m_source][(*ei).m_target].clear();
    }
    next_path_id = 0;

    // if(verbose)
    //     std::cout<<"step 1: cleared link load and path ids"<<std::endl;
    
    // clear and update routing table, update link load, update path ids
    for (int s = 0; s < num_vertices; s++) {
        for (int d = 0; d < num_vertices; d++) {
            routing_table[s][d].clear();
            if(s==d){
                continue;
            }
            std::vector<std::vector<Vertex>> paths = all_paths_all_s_d[s][d];
            // if(verbose)
            //     std::cout<<"step 1: cleared routing table for s = "<<s<<" d = "<<d<<std::endl;

            float ECMP_weight = 1.0/paths.size();
            for (std::vector<Vertex> path: paths){
                path_id path_id = next_path_id++;
                path_id_to_path[path_id] = path;
                routing_table[s][d][path_id] = ECMP_weight;
                for (int l = 0; l < path.size() - 1; l++) {
                    Vertex u = path[l];
                    Vertex v = path[l+1];
                    Edge e = boost::edge(u, v, G).first; // the boost::edge(u, v, G) returns a pair<Edge edge(u,v), bool found>
                    link_path_ids[e.m_source][e.m_target].push_back(path_id);
                    link_load[e.m_source][e.m_target] += ECMP_weight*M_R[s][d]/Cap_core; // TODO: to further optimize, divide Cap_core on the M_R at the beginning?
                    boost::put(boost::edge_weight, G, e, _alpha + pow(link_load[e.m_source][e.m_target],_beta)); // TODO: alternatively, update the weights outside this loop, may lead to speedup
                }
            }
        }
    }

    // if(verbose)
    //     std::cout<<"step 1: updated routing table, update link load, path ids, and weights"<<std::endl;


    // find the max value of link_load (float** , of size num_vertices*num_vertices)
    float max_load = 0.0;
    // EdgeIterator ei, ei_end;
    for (tie(ei, ei_end) = boost::edges(G); ei != ei_end; ei++) {
        if (link_load[(*ei).m_source][(*ei).m_target] > max_load) {
            max_load = link_load[(*ei).m_source][(*ei).m_target];
        }
    }
    // if(verbose)
    //     std::cout<<"step 1: found max link load"<<std::endl;

    // float max_load = 0;
    // for (tie(ei, ei_end) = boost::edges(G); ei != ei_end; ei++) {
    //     if (link_load[*ei] > max_load) {
    //         max_load = link_load[*ei];
    //     }
    // }
    result_max_loads_step_1.push_back(max_load);
    if(verbose)
        std::cout<<"result from step1: " << result_max_loads_step_1.back() << std::endl;


    // #ifdef DEBUG
    // std::cout << "Nexullance_IT step 1 done in debug mode " << std::endl;
    // #endif

    return;
}

bool Nexullance_IT::step_2(float _alpha, float _beta, float step, float threshold, int min_attempts, int max_attempts) {
    
    auto rng = std::default_random_engine {};
    int attempts = 0;
    std::list<float> max_loads_hist;
    long last_decreased_path_id = -1;
    long last_increased_path_id = -1;

    while (attempts < max_attempts) {
        float max_load = 0.0;
        std::vector<std::pair<int, int>> max_load_indices;

        EdgeIterator ei, ei_end;
        for (tie(ei, ei_end) = boost::edges(G); ei != ei_end; ei++) {
            size_t i = (*ei).m_source;
            size_t j = (*ei).m_target;
            if (link_load[i][j] > max_load) {
                max_load = link_load[i][j];
                max_load_indices.clear();
                max_load_indices.push_back({i, j});
            } else if (link_load[i][j] == max_load) {
                max_load_indices.push_back({i, j});
            }
        }
     
        max_loads_hist.push_back(max_load);
        result_max_load_step_2=max_load;

        if(max_load <= max_access_load){
            if (verbose){
                std::cout<<"max core link load reaches max access link load, terminating "<<std::endl;
            }
            num_attempts_step_2 += attempts;
            return false;
        }

        if(verbose)
            std::cout<<"step2, " << "step_value= "<< step << ", it= " << attempts << ", max_load = " << max_load << ",estimated phi = " 
                << total_flow/std::max(max_access_load, result_max_load_step_2)/(num_vertices*EPR)<< std::endl;   

        if(attempts > min_attempts){
            float recent_average_load = std::accumulate(std::prev(max_loads_hist.end(), min_attempts/2), max_loads_hist.end(), 0.0f)/((float)(min_attempts/2));
            if  ((recent_average_load - max_load)<threshold){
                if (verbose){
                    std::cout<<"step 2: low progress, terminating for step = "<< step <<std::endl;
                    std::cout<<"step 2: found max link load" << max_load <<std::endl;
                }
                num_attempts_step_2 += attempts;
                return true;
            }
        }
        bool success_attempt = false;

        // iterate through max_load_indices
        for (auto item: max_load_indices) {
            int i = item.first;
            int j = item.second;
            #ifdef DEBUG
            assert(link_load[i][j] == max_load && "link load should be max load");
            #endif
            std::vector<path_id> path_ids = link_path_ids[i][j];

            std::multimap<float, path_id, std::greater<float>> sorted_path_ids;
            for (auto old_path_id: path_ids) {
                std::vector<Vertex> old_path = path_id_to_path[old_path_id];
                int src = old_path.front();
                int dst = old_path.back();
                float contribution = routing_table[src][dst][old_path_id]*M_R[src][dst]/Cap_core;
                sorted_path_ids.insert(std::make_pair(contribution, old_path_id));
            }


            if (verbose) {
                std::cout << "Max load edge: (" << i << ", " << j << "), "
                          << "number of paths through this edge: " << sorted_path_ids.size() << std::endl;
            }

            for (auto item: sorted_path_ids) {
                path_id old_path_id = item.second;
                std::vector<Vertex> old_path = path_id_to_path[old_path_id];
                Vertex src = old_path.front();
                Vertex dst = old_path.back();

                // if(verbose){
                //     std::cout<<"step 2: starting with old path: " ;
                //     for(auto v: old_path)
                //         std::cout<<v<<" ";
                //     std::cout<<std::endl;
                // }
                        
                std::vector<std::vector<Vertex>> all_paths;
                compute_all_shortest_paths_single_s_d(G, src, dst, all_paths, weightmap);

                // if(verbose)
                //     std::cout<<"step 2: found " << all_paths.size() << " new paths for s = "<<src<<" d = "<<dst<<std::endl;

                std::shuffle(std::begin(all_paths), std::end(all_paths), rng);

                for(std::vector<Vertex> new_path: all_paths){
                    if(new_path != old_path){

                        float new_path_max_load = 0.0;
                        for (size_t l = 0; l < new_path.size() - 1; l++)
                            if (link_load[new_path[l]][new_path[l+1]] > new_path_max_load)
                                new_path_max_load = link_load[new_path[l]][new_path[l+1]];
                        
                        if(new_path_max_load == max_load)  continue;
                        #ifdef DEBUG
                        assert(new_path_max_load < max_load && "new path max load should be less than max load");
                        #endif

                        // if(verbose){
                        //     std::cout<<"step 2: starting with new path: " ;
                        //     for(auto v: new_path)
                        //         std::cout<<v<<" ";
                        //     std::cout<<std::endl;
                        // }
                        
                        //update the paths    
                        float delta_weigth = -1.0;
                        std::unordered_map<path_id, float>& current_routing_table = routing_table[src][dst];

                        auto iter = current_routing_table.find(old_path_id);
                        #ifdef DEBUG
                        assert(iter != current_routing_table.end());
                        #endif
                        float old_path_weight = iter->second;

                        bool new_path_found = false;
                        long new_path_id = -1;
                        float prev_path_weight = -1.0;
                        for (auto it = current_routing_table.begin(); it != current_routing_table.end(); ++it) {
                            if (path_id_to_path[it->first] == new_path) {
                                new_path_found = true;
                                new_path_id = it->first;
                                prev_path_weight = it->second;
                                break;
                            }
                        }
                        if(new_path_found){
                            #ifdef DEBUG
                            assert(new_path_id != -1);
                            assert(prev_path_weight != -1.0);
                            #endif
                            if ((old_path_id == last_increased_path_id) && (new_path_id == last_decreased_path_id)){
                                if(verbose){
                                    std::cout<<"step 2: stopped with the new path for avoiding deadlock "<<std::endl;
                                }
                                continue;
                            }
                            delta_weigth = std::min(std::min(step, std::min(old_path_weight, 1 - prev_path_weight)), 
                                            Cap_core*(max_load-new_path_max_load)/M_R[src][dst]); 

                            current_routing_table.erase(new_path_id);
                            current_routing_table.insert(std::make_pair(new_path_id, prev_path_weight+delta_weigth));

                        }else{
                            #ifdef DEBUG
                            assert(prev_path_weight == -1.0);
                            #endif

                            new_path_id = next_path_id++;
                            path_id_to_path[new_path_id] = new_path;
                            delta_weigth = std::min(std::min(step, old_path_weight), 
                                            Cap_core*(max_load-new_path_max_load)/M_R[src][dst]);
                            current_routing_table.insert(std::make_pair(new_path_id, delta_weigth));                
                            // current_routing_table[new_path_id] = delta_weigth;
                        }
                        success_attempt = true;
                        attempts++;
                        #ifdef DEBUG
                        assert(delta_weigth != -1.0);
                        #endif

                        last_decreased_path_id = old_path_id;
                        last_increased_path_id = new_path_id;

                        // iterate over the new path and update the link load and path ids
                        for (int l = 0; l < new_path.size() - 1; l++) {
                            Vertex u = new_path[l];
                            Vertex v = new_path[l+1];
                            Edge e = boost::edge(u, v, G).first; // the boost::edge(u, v, G) returns a pair<Edge edge(u,v), bool found>
                            link_load[e.m_source][e.m_target] += delta_weigth*M_R[src][dst]/Cap_core; // TODO: to further optimize, divide Cap_core on the M_R at the beginning?
                           
                            boost::put(boost::edge_weight, G, e, _alpha + pow(link_load[e.m_source][e.m_target] ,_beta));
                            if(!new_path_found)
                                link_path_ids[e.m_source][e.m_target].push_back(new_path_id);
                        }
                        current_routing_table.erase(old_path_id);
                        current_routing_table.insert(std::make_pair(old_path_id, old_path_weight-delta_weigth));

                        // iterate over the old path and update the link load and path ids
                        for (int l = 0; l < old_path.size() - 1; l++) {
                            Vertex u = old_path[l];
                            Vertex v = old_path[l+1];
                            Edge e = boost::edge(u, v, G).first; // the boost::edge(u, v, G) returns a pair<Edge edge(u,v), bool found>
                            link_load[e.m_source][e.m_target] -= delta_weigth*M_R[src][dst]/Cap_core; // TODO: to further optimize, divide Cap_core on the M_R at the beginning?
                            boost::put(boost::edge_weight, G, e, _alpha + pow(link_load[e.m_source][e.m_target],_beta));
                        
                            if (current_routing_table[old_path_id] < 0.00001) {
                                auto iter = std::find(link_path_ids[e.m_source][e.m_target].begin(), link_path_ids[e.m_source][e.m_target].end(), old_path_id);
                                if(iter != link_path_ids[e.m_source][e.m_target].end()){
                                    link_path_ids[e.m_source][e.m_target].erase(iter);
                                }
                            }
                        }
                        if (current_routing_table[old_path_id] < 0.00001) {
                            current_routing_table.erase(old_path_id);
                            path_id_to_path.erase(old_path_id);
                        }
                            
                        break;
                    }
                }

                if(success_attempt ){
                    break;
                }else{
                    continue;
                }

            }

            if(success_attempt){
                break;
            }else{
                continue;
            }
        }

        if(!success_attempt){
            if(verbose){
                std::cout<<"step 2: found max link load" << max_load <<std::endl;
                std::cout<<"step 2: no progress, terminating after " << attempts << " attempts"<<std::endl;
                }
            result_max_load_step_2=max_load;
            num_attempts_step_2 += attempts;
            return false;
        }
    }

    // if(verbose)
        // std::cout<<"step 2: found max link load" << max_load <<std::endl;
    std::cout<<"step 2: max number of attemtps reached with step = "<<step<<", threshold = "<<threshold<<", min_attempts = "<<min_attempts<<", max_attempts = "<<max_attempts<<std::endl;
    num_attempts_step_2 += attempts;
    return true;
}

void Nexullance_IT::optimize(int num_step_1, float alpha_step_1, float beta_step_1, int max_num_step_2, 
        float alpha_step_2, float beta_step_2, int step_2_min_attempts, float step_2_stepping_threshold, int step_2_max_attempts){
    
    assert(num_step_1 >= 1 and "num_step_1 should be a positive integer greater than 1.");
    if (verbose)    std::cout<<"Nexullance_IT starts optimization for the input demand matrix"<<std::endl; 
    
    for (int i = 0; i < num_step_1; i++) {
        step_1(alpha_step_1, beta_step_1);
    }
    float step = 0.5;
    for (int i = 0; i < max_num_step_2; i++) {
        if(step_2(alpha_step_2, beta_step_2, step, step_2_stepping_threshold, step_2_min_attempts, step_2_max_attempts)){
            step *= 0.5;
            continue;
        }else{
            break;
        }
    }
}


result_routing_table Nexullance_IT::get_routing_table(){
    result_routing_table result = result_routing_table();
    // iterate over "routing_table" and convert it to "result_routing_table"
    for (int s = 0; s < num_vertices; s++) {
        for (int d = 0; d < num_vertices; d++) {
            if (s==d)
                continue;

            std::vector< std::pair< std::vector<Vertex>,float> > paths = std::vector< std::pair< std::vector<Vertex>,float> >();
            for (auto item: routing_table[s][d]) {
                std::vector<Vertex> path = path_id_to_path[item.first];
                float weight = item.second;
                paths.push_back(std::make_pair(path, weight));
            }
            result.insert(std::make_pair(std::make_pair(s,d), paths));
        }
    }
    return result;
}

// float Nexullance_IT::get_average_path_length(){
//     float ave = 0;
//     for (int s = 0; s < num_vertices; s++) {
//         for (int d = 0; d < num_vertices; d++) {
//             if (s==d)
//                 continue;
//             for (auto item: routing_table[s][d]) {
//                 std::vector<Vertex> path = path_id_to_path[item.first];
//                 float weight = item.second;
//                 ave += path.size()*weight;
//             }
//         }
//     }
//     ave = ave/num_vertices/num_vertices;
//     return ave;
// }

float Nexullance_IT::get_max_core_load(){
    return result_max_load_step_2;
}

float Nexullance_IT::get_phi(){
    return total_flow/std::max(max_access_load, result_max_load_step_2)/(num_vertices*EPR);
}
