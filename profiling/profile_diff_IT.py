import sys

import topoResearch.global_helpers as gl
from topoResearch.topo_paths import REPO_ROOT
import numpy as np
import csv
import argparse
from topoResearch.nexullance.ultility import nexullance_exp_container
from diff_nexu.sampled_demand_analyzer import sampled_demand_analyzer

Cap_core = 10 #GBps
Cap_access = 10 #GBps

benchmark = ("FFT3D", 256, lambda size: f" nx={size} ny={size} nz={size} npRow=12")
    # ("Alltoall", 64, lambda size: f" bytes={size}")
    # ("Allreduce", 2048, lambda size: f" iterations=10 count={size}")

def main():
    parser = argparse.ArgumentParser(description='Profile IT with configurable number of samples')
    parser.add_argument('--num-samples', type=int, default=1, help='Number of samples to process (default: 1)')
    args = parser.parse_args()
    
    _num_samples = args.num_samples

    config = gl.ddf_configs[0]
    V = config[0]
    D = config[1]
    EPR = (D+1)//2
    topo_name="DDF"
    sampling_method = "sent"
    sampling_method = "enroute"

    bench, problem_size, args_func = benchmark
    args_str = args_func(problem_size)
    pickle_file = f"{REPO_ROOT}/diff_nexu/data/{bench}{args_str}_({V},{D})RRG_ECMP_ASP_{sampling_method}.pickle.gz"
    demand_analyzer = sampled_demand_analyzer(pickle_file)
    
    print(f"Running experiments for {_num_samples} samples...")
    
    M_EPs_s = demand_analyzer.get_inter_EP_demand_matrices(num_samples=_num_samples)[0]

    # now start with MD_Nexullance:
    exp_container = nexullance_exp_container(topo_name, V, D, EPR)
    exp_container.run_diff_nexullance_IT_batch_mode(M_EPs_s, auto_scaling=True)

    print(f"Completed experiments for {_num_samples} samples")


if __name__ == '__main__':
    main()