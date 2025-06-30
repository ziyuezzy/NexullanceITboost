#ifndef DIFF_IT_HPP
#define DIFF_IT_HPP

#include "definitions.hpp"
#include <boost/unordered/unordered_map.hpp>
// #include <list>
#include <unordered_map>
#include <tuple>

using path_id = size_t;

class diff_Nexullance_IT{

    public:
        diff_Nexullance_IT(Graph& _input_graph, const int _EPR, const float _Cap_core, 
            const float _Cap_access, const bool _verbose=false);
        ~diff_Nexullance_IT();

        // the idea of this class is to implement the diff-Nexullance algorithm, which is an on-line algorithm.
        // When exported to python, the instance of this class should be able to accept a new demand matrix at any time,
        // and optimize the routing table accordingly.

        IT_outputs optimize_for_M_EPs(float** M_EPs, float _alpha, float _beta, float threshold, 
            size_t max_num_step, int min_attempts, int max_attempts);

        // return <bool continue or not, num_attempts, max_core_load>
        std::tuple<bool, size_t, float> optimize_for_M_R_fixed_step(float** M_R, float _alpha, 
            float _beta, float threshold, float step, int min_attempts, int max_attempts, float max_access_load, float total_flow);

        // The class also accepts a list of demand matrices (for testing and analyzing), 
        // and optimize the routing table for each of them as if they are given consecutively.
        // void run_for_M_EPs_s(std::vector<float**> M_EPs_s, float _alpha, float _beta, float threshold, 
        //     float min_step, int min_attempts, int max_attempts);

        void initialize_routing_table();

    private:
        Graph G;
        const int EPR;
        size_t num_edges;
        size_t num_vertices;
        size_t M; // number of demand matrices as input
        
        const float Cap_core;
        const float Cap_access;
        const bool verbose;

        // void cleanup_for_next_matrix();
        // variable to re-used for each demand matrix
        std::vector<std::vector<Vertex>>** all_paths_all_s_d;
        property_map< Graph, edge_weight_t >::type weightmap;
        path_id next_path_id;
        std::unordered_map<path_id, std::vector<Vertex>> path_id_to_path; // using vector here, because the shortest-path algorithm return vector<Vertex> as a path
        std::unordered_map<path_id, float>** routing_table; // a 2D-array of map, first index corresponds to source switch id, second index the destination switch id.

        std::list<std::unordered_map<path_id, float>**> routing_tables_bin;

        result_routing_table export_routing_table(std::unordered_map<path_id, float>** table){
            result_routing_table result = result_routing_table();
            // iterate over "routing_table" and convert it to "result_routing_table"
            for (int s = 0; s < num_vertices; s++) {
                for (int d = 0; d < num_vertices; d++) {
                    if (s==d)
                        continue;

                    std::vector< std::pair< std::vector<Vertex>,float> > paths = std::vector< std::pair< std::vector<Vertex>,float> >();
                    for (auto item: table[s][d]) {
                        std::vector<Vertex> path = path_id_to_path[item.first];
                        float weight = item.second;
                        paths.push_back(std::make_pair(path, weight));
                    }
                    result.insert(std::make_pair(std::make_pair(s,d), paths));
                }
            }
            return result;
        }

        // std::unordered_map<path_id, float>** deep_copy_routing_table(std::unordered_map<path_id, float>** input_routing_table) {
        //     std::unordered_map<path_id, float>** copy = new std::unordered_map<path_id, float>*[num_vertices];
        //     for (size_t i = 0; i < num_vertices; i++) {
        //         copy[i] = new std::unordered_map<path_id, float>[num_vertices];
        //         for (size_t j = 0; j < num_vertices; j++)
        //             copy[i][j] = input_routing_table[i][j]; // std::unordered_map has a built-in copy constructor
        //     }
        //     routing_tables_bin.push_back(copy);
        //     return copy;
        // }

        float** link_load;
        // a mapping from link to path_ids
        std::vector<path_id>** link_path_ids;

        std::vector<IT_outputs> IT_outputs_list;
};
#endif
