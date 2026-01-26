import sys

from topoResearch.paths import IT_boost_debug, IT_boost_profile
from topoResearch.topo_paths import REPO_ROOT
sys.path.append(IT_boost_profile)
from Nexullance_IT_cpp import Nexullance_IT_fast_interface

import topoResearch.global_helpers as gl
import topoResearch.topologies.RRG as RRG
import numpy as np
import networkx as nx

import csv

configs = [(36, 5), (49, 6), (64, 7), (81, 8), (100, 9), (121, 10), (144, 11), (169, 12)]
# configs = [(36, 5), (49, 6)]

Cap_core = 10 #GBps
Cap_access = 10 #GBps

def main():


    # Create CSV file for results
    csv_filename = f"profiling_results_IT_fast_init.csv"
    with open(csv_filename, 'w', newline='') as csvfile:
        csvwriter = csv.writer(csvfile)
        # Write header
        csvwriter.writerow(['Nodes (V)', 'Degree (D)', 'Init Time (s)', 'RAM Usage (MB)'])
        
        # Process each network configuration
        for config in configs:
            V = config[0]
            D = config[1]
            _network = RRG.RRGtopo(V, D)
            arcs = _network.generate_graph_arcs()

            nexu_it_fast = Nexullance_IT_fast_interface(V, arcs, 10.0, 10.0, True)
            time, ram = nexu_it_fast.profile_initialization()
            ram_MB = ram/1024.0/1024.0
            
            # Write results to CSV
            csvwriter.writerow([V, D, time, ram_MB])
            csvfile.flush()  # write immediately

            # Print progress
            print(f"Processed network with {V} nodes and degree {D}: {time} seconds, {ram_MB:.2f} MB")

    # # use networkx to validate the results
    # _network.pre_calculate_APST_n(4)
    # APST_4_dict = _network.__getattribute__(f"APST_{4}")
    # # construct a validating routing table:
    # validating_RT = {}
    # for src in range(V):
    #     for dst in range(V):
    #         if src == dst:
    #             continue
    #         paths = APST_4_dict[(src, dst)] # contains multiple paths
    #         # Find the shortest path length
    #         min_path_length = min(len(path) for path in paths)
    #         # Count how many paths are of the shortest length
    #         shortest_paths = [path for path in paths if len(path) == min_path_length]
    #         num_shortest_paths = len(shortest_paths)
    #         # Compute the weight for shortest paths
    #         shortest_path_weight = 1.0 / num_shortest_paths
    #         validating_RT[(src, dst)] = []
    #         for path in paths:
    #             # Assign weight only to shortest paths
    #             weight = shortest_path_weight if len(path) == min_path_length else 0.0
    #             validating_RT[(src, dst)].append((list(path), weight))

    # # Compare the result with the validating routing table
    # for src in range(V):
    #     for dst in range(V):
    #         if src == dst:
    #             continue
    #         if (src, dst) not in result_RT:
    #             print(f"Missing route from {src} to {dst} in the result")
    #             continue
    #         result_paths = result_RT[(src, dst)]
    #         validating_paths = validating_RT[(src, dst)]
    #         if len(result_paths) != len(validating_paths):
    #             print(f"Different number of paths from {src} to {dst}: "
    #                   f"{len(result_paths)} vs {len(validating_paths)}")
    #             continue
    #         for i, (result_path, result_weight) in enumerate(result_paths):
    #             validating_path, validating_weight = validating_paths[i]
    #             if result_path != validating_path or round(result_weight, 3) != round(validating_weight, 3):
    #                 print(f"Path mismatch for {src} to {dst}: "
    #                       f"Result: {result_path}, Weight: {result_weight} "
    #                       f"Validating: {validating_path}, Weight: {validating_weight}")
    # print("Validation complete.")

if __name__ == '__main__':
    main()