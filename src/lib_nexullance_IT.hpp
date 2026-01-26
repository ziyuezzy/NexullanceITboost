#ifndef LIB_NEXULLANCE_IT_HPP
#define LIB_NEXULLANCE_IT_HPP

#include "definitions.hpp"
#include "graph_utility.hpp"
#include "Nexullance_IT.hpp"
#include "MD_Nexullance_IT.hpp"
#include <Eigen/Dense>

// Define a struct that contains all outputs from a run of Nexullance_IT, including the profiled time (RAM measurement will be in python calls), 
// the resulting max link load, the resulting phi and the routing table
struct IT_outputs{
    IT_outputs(double _elapsed_time, float _max_link_load, float _phi, result_routing_table _routing_table) : 
            elapsed_time(_elapsed_time), max_link_load(_max_link_load), phi(_phi), routing_table(_routing_table) { }
    double get_elapsed_time() const {return elapsed_time;}
    float get_max_core_link_load() const {return max_link_load;}
    float get_phi() const {return phi;}
    result_routing_table get_routing_table() const {return routing_table;}

    double elapsed_time;
    float max_link_load;
    float phi;
    result_routing_table routing_table;
};

// Similarly, define a struct for MD_Nexullance_IT. This includes the profiled time, 
// the resulting list of max link load, the list of phi and the routing table.
struct MD_IT_outputs{
    MD_IT_outputs(double _elapsed_time, std::vector<float> _max_link_loads, std::vector<float> _phis, result_routing_table _routing_table, float _Objective_func) : 
            elapsed_time(_elapsed_time), max_link_loads(_max_link_loads), phis(_phis), routing_table(_routing_table), Objective_func(_Objective_func) { }
    double get_elapsed_time() const {return elapsed_time;}
    std::vector<float> get_max_core_link_loads() const {return max_link_loads;}
    std::vector<float> get_phis() const {return phis;}
    result_routing_table get_routing_table() const {return routing_table;}
    float get_obj() const {return Objective_func;}

    double elapsed_time;
    std::vector<float> max_link_loads;
    std::vector<float> phis;
    result_routing_table routing_table;
    float Objective_func;
};

// TODO: a path can be simlified as "Vertex*"

class Nexullance_IT_interface {
    public:
    Nexullance_IT_interface(int V, Eigen::MatrixX2i arcs, const float Cap_core = 10 /*Gbps*/,
                                const float Cap_access = 10 /*Gbps*/, bool debug = false);
    ~Nexullance_IT_interface();

    Graph G;
    bool _debug = false;

    IT_outputs run_IT(Eigen::MatrixXf M_EPs, const int EPR);
    MD_IT_outputs run_MD_IT(std::vector<Eigen::MatrixXf> M_EPs_s, std::vector<float> weights, const int EPR);

    inline void set_parameters(float alpha, float beta, float stepping_threshold, 
                       int num_steppings, int max_attempts, int min_attempts, bool cal_least_margin=false){
        _alpha = alpha;
        _beta = beta;
        // _init_step = init_step;
        _stepping_threshold = stepping_threshold;
        _num_steppings = num_steppings;
        _max_attempts = max_attempts;
        _min_attempts = min_attempts;
        _cal_least_margin = cal_least_margin;
        }


    private:
    int _V;
    const float _Cap_core;
    const float _Cap_access;

    // For now, assume 'step 1' is only done once, with alpha=beta=1.0
    // So the only 'step 2' parameters are needed
    float _alpha=0.1f;
    float _beta=7.0f;
    // float _init_step = 0.5f;
    float _stepping_threshold = 0.0001f;
    int _num_steppings = 5; // number of times to adjust step size
    int _max_attempts = 1000000;
    int _min_attempts = 100;
    bool _cal_least_margin = false;

    // These above parameters are configured in method "set_parameters"
    
};


#endif