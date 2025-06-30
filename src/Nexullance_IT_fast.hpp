#ifndef NEXULLANCE_IT_FAST_HPP
#define NEXULLANCE_IT_FAST_HPP

#include "definitions.hpp"
#include <boost/unordered/unordered_map.hpp>
// #include <list>
#include <unordered_map>

// abbreviations for clarity of the code:
// using path_id = size_t;
typedef std::list<Vertex> path;
typedef std::pair<const path, float> path_info; // pair of the path and its weight

class Nexullance_IT_fast{
    public:
        Nexullance_IT_fast(Graph& _input_graph, const float _Cap_core,
                           const float _Cap_access, const bool _verbose=false);
        ~Nexullance_IT_fast();

        // find a paths shorter or equal to a certain length, construct mappings.
        // As the result of this method, 'routing_table', 'Map_passby', and 'Map_NOT_passby' will be filled correctly.
        // Input args:
        // 1. the first argument is the maximum path length that will be found in the graph and put into the routing table
        // 2. the second argument sets limitation to setting weights at this initalization stage: 
            // only paths shorter or equal to 'initial_weight_max_path_length' will be assigned equal weights 
            // (one divided by the number of paths that are shorter or equal to 'initial_weight_max_path_length' between this s-d pair).
            // for all paths that are longer than this value, the weight will be set to zero.
        void initialization(size_t max_path_length, size_t initial_weight_max_path_length);

        IT_outputs optimize_for_M_EPs(float** M_EPs, float threshold, 
            size_t max_num_step, int min_attempts, int max_attempts);

        // return <bool continue or not, num_attempts, max_core_load>
        std::tuple<bool, size_t, float> optimize_for_M_R_fixed_step(float** M_R, float threshold, 
            float step, int min_attempts, int max_attempts, float max_access_load, float total_flow);

        result_routing_table get_routing_table();
        ulong get_RAM_after_init();
    private:
        Graph G;
        size_t num_edges;
        size_t num_vertices;
        
        const float Cap_core;
        const float Cap_access;
        const bool verbose;

        // some data structure to store paths and their mappings

        // "routing_table" is a 3D vector, first index corresponds to source switch id, 
        // second index the destination switch id, third dimension of the vector is a list of paths between the source and destination,
        // each element is a pointer to a path_info (a pair of the path and its weight).
        std::vector< std::vector< std::vector< path_info* > > > routing_table;

        // a mapping from link (edge in the graph) to paths that pass by this link.
        boost::unordered_map<Edge, std::vector<path_info*>> Map_passby; 
        // Alternatively, the vector can be a ordered_map of path sorted by their load contributions (demand*weight)
        // However, then the ordered map need to be updated whenever the path weights are changed, which might be more expensive than keeping an unsorted vector?

        // a mapping from link to paths that DO NOT pass by this link.
        boost::unordered_map<Edge, std::vector< std::vector< std::vector< path_info*> > > > Map_NOT_passby;
        // The value of the map is a 3-D vector, where the first index corresponds to the source switch id,
        // the second index the destination switch id, and the third dimension is a vector of path_info pointers that do not pass by this link.
        // This is used to quickly find paths that do not pass by a certain link.
        // TODO: 2. it can also be std::vector< std::vector<  std::unordered_map<Edge, std::vector<path_info*> > > >, which is better?

};

#endif
