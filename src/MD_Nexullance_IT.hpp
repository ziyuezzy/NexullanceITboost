#ifndef MD_NEXULLANCE_IT_HPP
#define MD_NEXULLANCE_IT_HPP

#include "definitions.hpp"
#include <boost/unordered/unordered_map.hpp>
// #include <list>
#include <unordered_map>
#include <map>

using path_id = size_t;

class MD_Nexullance_IT{
    // _nx_graph: nxGraph, _M_R: np.ndarray, _Cap_core: float, _verbose:bool=False
    public:
        MD_Nexullance_IT(Graph& _input_graph, std::vector<float**> _M_EPs_s, const std::vector<float> _M_weights, 
                        const int _EPR, const float _Cap_core, const float _Cap_access, const bool _verbose=false);
        ~MD_Nexullance_IT();

        void get_finial_max_load();
        void step_1(float _alpha, float _beta);
        bool step_2(float _alpha, float _beta, float step, float threshold=0.001, int min_attempts=50, int max_attempts=100000, bool cal_least_margin = false);
        void optimize(int num_step_1, float alpha_step_1, float beta_step_1, int max_num_step_2, float alpha_step_2, float beta_step_2, 
                        int method_2_min_attempts, float method_2_threshold, int method_2_max_attempts, bool cal_least_margin);

        // std::list<float> result_weighted_phi_step_1;
        // std::list<float> result_max_loads_step_2;
        size_t num_attempts_step_2 = 0;

        result_routing_table get_routing_table();
        // float get_average_path_length();
        inline std::vector<float> get_max_core_load_vec(){
            return max_core_loads;
        }
        inline std::vector<float> get_phis(){
            return phis;
        }
        // inline float get_weighted_sum_phi(){
        //     return weighted_sum_phi;
        // }
        inline float get_obj(){
            return 1/Objective_func_inverse;
        }

    private:
        Graph G;
        const int EPR;
        size_t num_edges;
        size_t num_vertices;
        size_t M; // number of demand matrices as input
        
        std::vector<float**> M_Rs;
        const std::vector<float> M_weights;
        std::vector<float> cumulative_weights;
        const float Cap_core;
        const float Cap_access;
        const bool verbose;
        
        std::vector<float> max_access_loads; // max access link load under each demand matrix
        std::vector<float> total_flows_per_EP; // total flow per EP in the demand matrix, under each demand matrix
        std::vector<std::vector<std::pair<size_t, size_t>>> max_core_load_links; // to store the max-loaded core link for each demand matrix

        std::vector<float> max_core_loads; // max core load under each demand matrix
        std::vector<float> phis;
        // float weighted_sum_phi;
        float Objective_func_inverse;


        path_id next_path_id;
        std::unordered_map<path_id, std::vector<Vertex>> path_id_to_path; // using vector here, because the shortest-path algorithm return vector<Vertex> as a path
        std::unordered_map<path_id, float>** routing_tables; // a 2D-array of map, first index corresponds to source router id, second index the destination router id.
        std::vector<boost::unordered_map<Edge, float>> link_load_vec; // link load for each demand matrix
        boost::unordered_map<Edge, std::vector<path_id>> link_path_ids;
        // float final_max_load;
        
        std::vector<std::vector<Vertex>>** all_paths_all_s_d;
        property_map< Graph, edge_weight_t >::type weightmap;

};
#endif // MD_NEXULLANCE_IT_HPP
