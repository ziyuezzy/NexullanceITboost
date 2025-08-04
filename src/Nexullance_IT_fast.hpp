#ifndef NEXULLANCE_IT_FAST_HPP
#define NEXULLANCE_IT_FAST_HPP

#include "definitions.hpp"
#include "IT_fast_definitions.hpp"
#include <boost/unordered/unordered_map.hpp>
// #include <list>
#include <unordered_map>
#include <random>

// // abbreviations for clarity of the code: (moved to IT_fast_definitions.hpp)
// typedef std::list<Vertex> path;
// typedef std::pair<const path, float> path_info; // pair of the path and its weight

class Nexullance_IT_fast{
    public:
        Nexullance_IT_fast(Graph& _input_graph, const float _Cap_core,
                           const float _Cap_access, const bool _verbose=false);
        ~Nexullance_IT_fast();

        // find a paths shorter or equal to a certain length, construct mappings.
        // As the result of this method, 'initial_RT', 'Map_passby', and 'Map_NOT_passby' will be filled correctly.
        // Input args:
        // 1. the first argument is the maximum path length that will be found in the graph and put into the routing table
        // 2. the second argument sets limitation to setting weights at this initalization stage: 
            // only paths shorter or equal to 'initial_weight_max_path_length' will be assigned equal weights 
            // (one divided by the number of paths that are shorter or equal to 'initial_weight_max_path_length' between this s-d pair).
            // for all paths that are longer than this value, the weight will be set to zero.
        void initialization(size_t max_path_length, size_t initial_weight_max_path_length, bool disjoint_paths = false);

        IT_outputs optimize_for_M_EPs(float** M_EPs, size_t EPR, size_t max_num_fix_step, bool from_initial_RT);

        // return <bool continue or not, num_attempts, max_core_load>
        std::tuple<bool, size_t, float> optimize_for_M_R_fixed_step(float** M_R, float step, float max_access_load, float total_flow);

        ulong get_RAM_after_init();

        void set_params(size_t _min_attempts, size_t _max_attempts, float _progress_threshold) {
            min_attempts = _min_attempts;
            max_attempts = _max_attempts;
            progress_threshold = _progress_threshold;
        }


        result_routing_table get_initial_routing_table(){
            return get_routing_table(initial_RT);
        }

        result_routing_table get_previous_result_routing_table(){
            return get_routing_table(previous_result_RT);
        }

    private:
        Graph G;
        size_t num_edges;
        size_t num_vertices;

        size_t min_attempts = 100; // minimum number of attempts (iterations) per fixed step
        size_t max_attempts = 1000000; // maximum number of attempts (iterations) per fixed step
        float progress_threshold = 0.0001f; // threshold for stopping a fixed step
        
        const float Cap_core;
        const float Cap_access;
        const bool verbose;

        result_routing_table get_routing_table(std::vector< std::vector< std::vector< path_info* > > > routing_table);

        // return the selected alternative path and the maximum load along the path
        std::pair<path_info*, float> select_alternative_path_random(path_info* old_path_info, Edge max_load_edge, float old_path_max_load);

        // TODO: another strategy for selecting alternative paths:
        // iterative over the paths that do not pass through the edge with max load, as long as the max load on the new path is lower than a value, we select it.
        // maybe this can be more efficient than the random selection strategy?

        // TODO: implement a selecting strategy that calculates the max load of all possible paths and selects the one with the lowest max load.
        // Use this to compare with the Dijkstra search?
        

        // TODO: implement Dijkstra search, but then the framework needs to be changed to support it.

        // some data structure to store paths and their mappings

        // a "routing table" is a 3D vector, first index corresponds to source switch id, 
        // second index the destination switch id, third dimension of the vector is a list of paths between the source and destination,
        // each element is a pointer to a path_info (a pair of the path and its weight).
        std::vector< std::vector< std::vector< path_info* > > > previous_result_RT;
        std::vector< std::vector< std::vector< path_info* > > > initial_RT;
        std::vector< std::vector< std::vector< path_info* > > > current_RT;

        // a mapping from link (edge in the graph) to paths that pass by this link, the paths are sorted by their flow (contribution to the load).
        boost::unordered_map<Edge, descend_loadcontri_pathinfos> Map_passby;
        // Alternatively, the vector can be a ordered_map of path sorted by their load contributions (demand*weight)
        // However, then the ordered map need to be updated whenever the path weights are changed, which might be more expensive than keeping an unsorted vector?

        // a mapping from link to paths that DO NOT pass by this link.
        boost::unordered_map<Edge, std::vector< std::vector< std::vector< path_info*> > > > Map_NOT_passby;
        // The value of the map is a 3-D vector, where the first index corresponds to the source switch id,
        // the second index the destination switch id, and the third dimension is a vector of path_info pointers that do not pass by this link.
        // This is used to quickly find paths that do not pass by a certain link.
        // TODO: 2. it can also be std::vector< std::vector<  std::unordered_map<Edge, std::vector<path_info*> > > >, which is better?

        // a mapping from a link to its load
        boost::unordered_map<Edge, float> link_load;

            
        std::default_random_engine rng;

};

#endif
