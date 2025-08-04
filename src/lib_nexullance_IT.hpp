#ifndef LIB_NEXULLANCE_IT_HPP
#define LIB_NEXULLANCE_IT_HPP

#include "definitions.hpp"
#include "graph_utility.hpp"
#include "Nexullance_IT.hpp"
#include "MD_Nexullance_IT.hpp"
#include "diff_Nexullance_IT.hpp"
#include "Nexullance_IT_fast.hpp"
// #include <Eigen/Dense>

// // for differential nexullance_IT, the input is a list of demand matrices, 
// // the result contains is a list of elapsed time,  a list of max link loads, a list of phi's and a list of routing tables
// struct diff_IT_outputs{
//     diff_IT_outputs(vector<float> _elapsed_time_list, vector<float> _max_link_load_list, 
//                     vector<float> _phi_list, vector<result_routing_table> _routing_table_list) : 
//             elapsed_time_list(_elapsed_time_list), max_link_load_list(_max_link_load_list), 
//             phi_list(_phi_list), routing_table_list(_routing_table_list) {}
//     public:
//     vector<float> get_elapsed_time_list() const {return elapsed_time_list;}
//     vector<float> get_max_link_load_list() const {return max_link_load_list;}
//     vector<float> get_phi_list() const {return phi_list;}
//     vector<result_routing_table> get_routing_table_list() const {return routing_table_list;}

//     private:
//     vector<float> elapsed_time_list;
//     vector<float> max_link_load_list;
//     vector<float> phi_list;
//     vector<result_routing_table> routing_table_list;
// };

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
                       int max_num_step2, int max_attempts, int min_attempts, bool cal_least_margins=false){
        _alpha = alpha;
        _beta = beta;
        // _init_step = init_step;
        _stepping_threshold = stepping_threshold;
        _max_num_step2 = max_num_step2;
        _max_attempts = max_attempts;
        _min_attempts = min_attempts;
        _cal_least_margins = cal_least_margins;
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
    int _max_num_step2 = 5; // number of times to adjust step size
    int _max_attempts = 1000000;
    int _min_attempts = 100;
    bool _cal_least_margins = false;

    // These above parameters are configured in method "set_parameters"
    
};

class diff_Nexullance_IT_interface {
    public:
    diff_Nexullance_IT_interface(int V, Eigen::MatrixX2i arcs, const float Cap_core = 10 /*Gbps*/,
                                const float Cap_access = 10 /*Gbps*/, bool online_mode = false, bool debug = false);
    ~diff_Nexullance_IT_interface();

    Graph G;
    bool _debug = false;

    // input a list of demand matrices, as if they are consecutively sampled demand matrices from a certain application
    std::vector<IT_outputs> run_for_batch_matrices(std::vector<Eigen::MatrixXf> M_EPs_s, const int EPR);

    // initialize for online mode, the demand matrices are not known in advance
    // new demand matrices will be added one by one
    void _initialize(const int EPR);
    IT_outputs add_next_matrix(Eigen::MatrixXf M_EPs);

    inline void set_parameters(float alpha, float beta, float stepping_threshold, 
                       int max_num_step2, int max_attempts, int min_attempts){
        _alpha = alpha;
        _beta = beta;
        // _init_step = init_step;
        _stepping_threshold = stepping_threshold;
        _max_num_step2 = max_num_step2;
        _max_attempts = max_attempts;
        _min_attempts = min_attempts;
        }


    private:

    diff_Nexullance_IT* diff_nexu_it;
    int _V;
    const bool _online_mode;
    int num_EPs;
    const float _Cap_core;
    const float _Cap_access;
    
    float _alpha=0.1f;
    float _beta=7.0f;
    // float _init_step = 0.5f;
    float _stepping_threshold = 0.0001f;
    int _max_num_step2 = 5; // number of times to adjust step size
    int _max_attempts = 1000000;
    int _min_attempts = 100;

    // These above parameters are configured in method "set_parameters"
    
};

class Nexullance_IT_fast_interface {
    public:
    Nexullance_IT_fast_interface(int V, Eigen::MatrixX2i arcs, const float Cap_core = 10 /*Gbps*/,
                                const float Cap_access = 10 /*Gbps*/, bool debug = false);
    ~Nexullance_IT_fast_interface();

    Graph G;
    bool _debug = false;

    void _initialize(size_t max_path_length, size_t initial_weight_max_path_length, bool disjoint_paths);
    
    // profile the calculation time and memory usage of the initialization method
    std::pair<double, size_t> profile_initialization(size_t max_path_length, size_t initial_weight_max_path_length);

    result_routing_table get_initial_routing_table();
    result_routing_table get_previous_result_routing_table();

    IT_outputs add_next_matrix(Eigen::MatrixXf M_EPs, size_t max_num_fix_step, bool from_initial_RT);

    // input a list of demand matrices, as if they are consecutively sampled demand matrices from a certain application
    std::vector<IT_outputs> run_for_batch_matrices(std::vector<Eigen::MatrixXf> M_EPs_s, bool from_initial_RT);

    void set_parameters(size_t _min_attempts, size_t _max_attempts, float _progress_threshold) {
        nexu_it_fast->set_params(_min_attempts, _max_attempts, _progress_threshold);
    }


    private:

    Nexullance_IT_fast* nexu_it_fast;
    int _V;
    // int num_EPs;
    const float _Cap_core;
    const float _Cap_access;
    
    
};

#endif