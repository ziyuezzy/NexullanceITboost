#ifndef GRAPH_DEFINITIONS_HPP
#define GRAPH_DEFINITIONS_HPP
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>
#include <vector>

// #include <Eigen/Dense>
using namespace boost;

// Define a directed graph with integer vertex index and float edge weights
typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::directedS,
                               boost::property<boost::vertex_index_t, int>,
                               boost::property<boost::edge_weight_t, float>
                              > Graph;

// // Define a directed graph with integer vertex index, float edge weights, and integer edge indices
// typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::directedS,
//                                boost::property<boost::vertex_index_t, int>,
//                                boost::property<boost::edge_weight_t, float>,
//                                boost::property<boost::edge_index_t, int>
//                               > Graph;


typedef boost::graph_traits<Graph>::vertex_descriptor Vertex;
typedef boost::graph_traits<Graph>::edge_descriptor Edge;
typedef boost::graph_traits<Graph>::edge_iterator EdgeIterator;
typedef boost::graph_traits<Graph>::vertex_iterator VertexIterator;
typedef boost::graph_traits<Graph>::out_edge_iterator OutEdgeIterator;


typedef std::map<std::pair<Vertex,Vertex>, std::vector< std::pair<std::vector<Vertex>,float> > > result_routing_table;

// void test_shortest_paths(std::string input_graph_path); // implementaion in the bin

// Define a struct that contains all outputs from a run of Nexullance_IT, including the profiled time (RAM measurement will be in python calls), 
// the resulting max link load, the resulting phi and the routing table
struct IT_outputs{
    IT_outputs(double _elapsed_time, float _max_link_load, float _phi, result_routing_table _routing_table, size_t _num_attempts) : 
            elapsed_time(_elapsed_time), max_core_link_load(_max_link_load), phi(_phi), routing_table(_routing_table), num_attempts(_num_attempts) { }
    public:
    double get_elapsed_time() const {return elapsed_time;}
    float get_max_core_link_load() const {return max_core_link_load;}
    float get_phi() const {return phi;}
    result_routing_table get_routing_table() const {return routing_table;}
    size_t get_num_attempts() const {return num_attempts;}
    
    private:
    double elapsed_time;
    float max_core_link_load;
    float phi;
    result_routing_table routing_table;
    size_t num_attempts;
};

// Similarly, define a struct for MD_Nexullance_IT. The profiled time for execution is returned, 
// the result also contains a list of max link load, a list of phi and the routing table.
struct MD_IT_outputs{
    MD_IT_outputs(double _elapsed_time, std::vector<float> _max_link_loads, std::vector<float> _phis, result_routing_table _routing_table, float _Objective_func, size_t _num_attempts) : 
            elapsed_time(_elapsed_time), max_core_link_loads(_max_link_loads), phis(_phis), routing_table(_routing_table), Objective_func(_Objective_func), num_attempts(_num_attempts) { }
    public:
    double get_elapsed_time() const {return elapsed_time;}
    std::vector<float> get_max_core_link_loads() const {return max_core_link_loads;}
    std::vector<float> get_phis() const {return phis;}
    result_routing_table get_routing_table() const {return routing_table;}
    float get_obj() const {return Objective_func;}
    size_t get_num_attempts() const {return num_attempts;}

    private:
    double elapsed_time;
    std::vector<float> max_core_link_loads;
    std::vector<float> phis;
    result_routing_table routing_table;
    float Objective_func;
    size_t num_attempts;
};

#endif // GRAPH_DEFINITIONS_HPP
