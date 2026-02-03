#ifndef NEXULLANCE_IT_HPP
#define NEXULLANCE_IT_HPP

#include "definitions.hpp"
#include <boost/unordered/unordered_map.hpp>
// #include <list>
#include <unordered_map>

using path_id = size_t;

class Nexullance_IT{
    // _nx_graph: nxGraph, _M_R: np.ndarray, _Cap_core: float, _verbose:bool=False
    public:
        Nexullance_IT(Graph& _input_graph, const float** _M_EPs, const int _EPR, const float _Cap_core, 
                        const float _Cap_access, const bool _verbose=false);
        ~Nexullance_IT();

        // void get_finial_max_load();
        void step_1(float _alpha, float _beta);
        bool step_2(float _alpha, float _beta, float step, float threshold=0.001, 
                    int min_attempts=50, int max_attempts=100000);
        void optimize(int num_step_1, float alpha_step_1, float beta_step_1, int max_num_step_2, 
                    float alpha_step_2, float beta_step_2, int step_2_min_attempts, int step_2_stepping_threshold, int step_2_max_attemtps);

        std::list<float> result_max_loads_step_1;
        // std::list<float> result_max_loads_step_2;
        float result_max_load_step_2;
        size_t num_attempts_step_2 = 0;

        result_routing_table get_routing_table();
        // float get_average_path_length();
        float get_phi() const;
        float get_max_core_load() const;

    private:
        Graph G;
        const int EPR;
        size_t num_edges;
        size_t num_vertices;
        float** M_R;
        const float Cap_core;
        const float Cap_access;
        const bool verbose;
        float max_access_load;
        float total_flow = 0.0f;

        path_id next_path_id;
        std::unordered_map<path_id, std::vector<Vertex>> path_id_to_path; // using vector here, because the shortest-path algorithm return vector<Vertex> as a path
        std::unordered_map<path_id, float>** routing_tables; // a 2D-array of map, first index corresponds to source router id, second index the destination router id.
        boost::unordered_map<Edge, float> link_load;
        boost::unordered_map<Edge, std::vector<path_id>> link_path_ids;
        // float final_max_load;
        
        std::vector<std::vector<Vertex>>** all_paths_all_s_d;
        property_map< Graph, edge_weight_t >::type weightmap;

        // Helper method to extract max link load
        float compute_max_link_load() const;
        
        // Helper method to update edge weight based on load
        void update_edge_weight(Edge e, float alpha, float beta);
};
#endif // NEXULLANCE_IT_HPP
