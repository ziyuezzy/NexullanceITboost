#include "diff_Nexullance_IT.hpp"
#include "graph_utility.hpp"
#include <boost/graph/adjacency_list.hpp>
// #include <boost/property_map/property_map.hpp>
#include <iostream>
#include <algorithm>
#include <numeric>
#include <iterator>
#include <random>
#include <chrono>

diff_Nexullance_IT::diff_Nexullance_IT(Graph& _input_graph, const int _EPR, const float _Cap_core, 
    const float _Cap_access, const bool _verbose): 
    G(_input_graph), EPR(_EPR), Cap_core(_Cap_core), Cap_access(_Cap_access), verbose(_verbose){

    num_edges = boost::num_edges(G);
    num_vertices = boost::num_vertices(G);   

    // initialize data structures
    all_paths_all_s_d = new std::vector<std::vector<Vertex>>*[num_vertices];
    for (int i = 0; i < num_vertices; i++) {
        all_paths_all_s_d[i] = new std::vector<std::vector<Vertex>>[num_vertices];
    }
    weightmap = get(edge_weight, G);

    link_load = new float*[num_vertices];
    for (int i = 0; i < num_vertices; i++)      link_load[i] = new float[num_vertices]{0.0f};

    link_path_ids = new std::vector<path_id>*[num_vertices];
    for (int i = 0; i < num_vertices; i++)      link_path_ids[i] = new std::vector<path_id>[num_vertices];
    
    routing_table = new std::unordered_map<path_id, float>*[num_vertices];
    for (int i = 0; i < num_vertices; i++)      routing_table[i] = new std::unordered_map<path_id, float>[num_vertices];

    initialize_routing_table();
}

diff_Nexullance_IT::~diff_Nexullance_IT() {
    for (int i = 0; i < num_vertices; i++) {
        delete[] link_load[i];
        delete[] link_path_ids[i];
        delete[] all_paths_all_s_d[i];
        delete[] routing_table[i];
    }
    delete[] link_load;
    delete[] link_path_ids;
    delete[] all_paths_all_s_d;
    delete[] routing_table;

    // TODO
    
}

void diff_Nexullance_IT::initialize_routing_table() {
    compute_all_shortest_paths_all_s_d(G, all_paths_all_s_d, weightmap);
    if(verbose)
        std::cout<<"diff nexu initialization: computed all shortest paths all s d"<<std::endl;

    // then set all loads to 0.0, empty path ids, empty path_id_to_path
    path_id_to_path.clear();
    EdgeIterator ei, ei_end;
    for (tie(ei, ei_end) = boost::edges(G); ei!= ei_end; ei++)
        link_path_ids[(*ei).m_source][(*ei).m_target].clear();

    next_path_id = 0;
    
    // clear and update routing table, assign path ids
    for (int s = 0; s < num_vertices; s++) {
        for (int d = 0; d < num_vertices; d++) {
            routing_table[s][d].clear();
            if(s==d)    continue;
            std::vector<std::vector<Vertex>> paths = all_paths_all_s_d[s][d];

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
                }
            }
        }
    }
    if(verbose)     std::cout<<"diff nexu initialization: updated routing table and assigned path ids"<<std::endl;
    return;
}

// void diff_Nexullance_IT::run_for_M_EPs_s(std::vector<float**> M_EPs_s, float _alpha, float _beta, 
//         float threshold, float min_step, int min_attempts, int max_attempts){
//     #ifdef DEBUG
//     assert(_online_mode == false);
//     #endif
//     size_t num_M = M_EPs_s.size();
//     for (size_t i = 0; i < num_M; i++){
//         IT_outputs result = optimize_for_M_EPs(M_EPs_s[i], _alpha, _beta, threshold, min_step, min_attempts, max_attempts);
//         IT_outputs_list.push_back(result);
//     }
// }

IT_outputs diff_Nexullance_IT::optimize_for_M_EPs(float** M_EPs, float _alpha, float _beta, float threshold, 
            size_t max_num_step2, int min_attempts, int max_attempts){

    if (verbose)    std::cout<<"diff_Nexullance_IT starts optimization for the input demand matrix"<<std::endl; 
    // start timer
    auto start = std::chrono::high_resolution_clock::now();
    // preparation
    float** M_R;
    std::pair<float, float> temp_result=procress_M_EPs(const_cast<const float**>(M_EPs), num_vertices, EPR, &M_R);
    float max_access_load = temp_result.first/Cap_access;
    float total_flow = temp_result.second;

    // calculate link loads and therefore the weights on graph edges
    EdgeIterator ei, ei_end;
    for (tie(ei, ei_end) = boost::edges(G); ei!= ei_end; ei++)
        link_load[(*ei).m_source][(*ei).m_target] = 0;

    for (int s = 0; s < num_vertices; s++) {
        for (int d = 0; d < num_vertices; d++) {
            for (auto it = routing_table[s][d].begin(); it != routing_table[s][d].end(); ++it) {
                path_id path_id = it->first;
                float ECMP_weight = it->second;
                std::vector<Vertex> path = path_id_to_path[path_id];
                for (int l = 0; l < path.size() - 1; l++) {
                    Vertex u = path[l];
                    Vertex v = path[l+1];
                    Edge e = boost::edge(u, v, G).first; // the boost::edge(u, v, G) returns a pair<Edge edge(u,v), bool found>
                    link_load[e.m_source][e.m_target] += ECMP_weight*M_R[s][d]/Cap_core;
                    boost::put(boost::edge_weight, G, e, _alpha + pow(link_load[e.m_source][e.m_target],_beta));
                }
            }
        }
    }

    bool to_continue = true;
    size_t tot_attempts = 0;
    size_t num_attempts;
    float max_core_load;
    float step = 0.5;
    size_t num_step2s = 0;
    while ((max_num_step2 > num_step2s) && to_continue)
    {
        std::tie(to_continue, num_attempts, max_core_load) = optimize_for_M_R_fixed_step(M_R, 
            _alpha, _beta, step, threshold, min_attempts, max_attempts, max_access_load, total_flow);
        step *= 0.5;
        tot_attempts += num_attempts;
        num_step2s++;
    }

    // free the matrix
    delete[] M_R;

    // end timer
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    if(verbose)     std::cout<<"elapsed time = " << elapsed.count() << " s" << std::endl;
    float phi = total_flow/std::max(max_access_load, max_core_load)/(num_vertices*EPR);
    
    IT_outputs result = IT_outputs(elapsed.count(), max_core_load, phi, export_routing_table(routing_table), tot_attempts);
    return result;
}


std::tuple<bool, size_t, float> diff_Nexullance_IT::optimize_for_M_R_fixed_step(float** M_R, float _alpha, 
    float _beta, float step, float threshold, int min_attempts, int max_attempts, float max_access_load, float total_flow){
    // calculate the load based on the new demand matrix
    // fine max loaded link and re-route the paths
    // the routing_table, link_path_ids and path_id_to_path should already be calculated (
    // based on the previous demand matrix or in the constructor)

    // start the main loop
    auto rng = std::default_random_engine {};
    size_t attempts = 0;
    std::list<float> max_loads_hist;
    float result_max_load = 0.0;
    long last_decreased_path_id = -1;
    long last_increased_path_id = -1;
    
    while (attempts < max_attempts) {

        // 2.1 find the max value of link_load
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
        result_max_load=max_load;

        if(max_load <= max_access_load){
            if (verbose){
                std::cout<<"max core link load reaches max access link load, terminating "<<std::endl;
            }
            return std::make_tuple(false, attempts, result_max_load);
        }

        if(verbose)
            std::cout<< "step_value= "<< step << ", it="<< attempts << ", max_load = " << result_max_load << ",estimated phi = " 
                << total_flow/std::max(max_access_load, result_max_load)/(num_vertices*EPR)<< std::endl;   
    
        if((attempts > min_attempts) && (( std::accumulate(std::prev(max_loads_hist.end(), min_attempts/2), max_loads_hist.end(), 0.0f)/((float)min_attempts) - max_load)<threshold)){
            if (verbose){
                std::cout<<"diff nexu:: low progress, terminating for step = "<< step <<std::endl;
                std::cout<<"diff nexu:: found max link load" << result_max_load <<std::endl;
            }
            return std::make_tuple(true, attempts, result_max_load);
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

            for (auto item: sorted_path_ids) {
                path_id old_path_id = item.second;
                std::vector<Vertex> old_path = path_id_to_path[old_path_id];
                Vertex src = old_path.front();
                Vertex dst = old_path.back();

                // if(verbose){
                //     std::cout<<"diff nexu:: starting with old path: " ;
                //     for(auto v: old_path)
                //         std::cout<<v<<" ";
                //     std::cout<<std::endl;
                // }
                        
                std::vector<std::vector<Vertex>> all_paths;
                compute_all_shortest_paths_single_s_d(G, src, dst, all_paths, weightmap);

                // if(verbose)
                //     std::cout<<"diff nexu:: found " << all_paths.size() << " new paths for s = "<<src<<" d = "<<dst<<std::endl;

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
                        //     std::cout<<"diff nexu:: starting with new path: " ;
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
                                    std::cout<<"diff nexu:: stopped with the new path for avoiding deadlock "<<std::endl;
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
                std::cout<<"found max link load" << max_load <<std::endl;
                std::cout<<"no progress, terminating after " << attempts << " attempts with step_value = "<< step <<std::endl;
                }
            return std::make_tuple(false, attempts, result_max_load);
        }
    }
         
    // if(verbose)
        // std::cout<<"diff nexu:: found max link load" << max_load <<std::endl;
    std::cout<<"max number of attemtps reached with step_value = "<<step<<", threshold = "<<threshold<<", min_attempts = "<<min_attempts<<", max_attempts = "<<max_attempts<<std::endl;
    return std::make_tuple(true, attempts, result_max_load);    
}