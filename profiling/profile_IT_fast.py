import sys
sys.path.append("/groups/ilabt-imec-be/hpcnetworksimulation/ziyzhang/EFM_experiments/")
from topoResearch.paths import IT_boost_debug, IT_boost_release
# sys.path.append(IT_boost_release)
sys.path.append(IT_boost_debug)

from Nexullance_IT_cpp import Nexullance_IT_fast_interface
import topoResearch.global_helpers as gl
import topoResearch.topologies.RRG as RRG
import numpy as np

V = 36
D = 5
EPR = (D+1)//2
_network = RRG.RRGtopo(V, D)
ASP, _ = _network.calculate_all_shortest_paths()
ECMP_ASP = gl.ECMP(ASP)
arcs = _network.generate_graph_arcs()

Cap_remote = 10 #GBps
Cap_local = 10 #GBps


print("shift-1 demand matrix:")

M_EPs = gl.generate_shift_traffic_demand_matrix(V, EPR, 1)
remote_link_flows, local_link_flows = _network.distribute_M_EPs_on_weighted_paths(ECMP_ASP, EPR, M_EPs)
max_remote_link_load = np.max(remote_link_flows)/Cap_remote
max_local_link_load = np.max(local_link_flows)/Cap_local
# adapt the traffic scaling factor to 10x saturation
traffic_scaling = 10.0/max(max_local_link_load, max_remote_link_load)
M_EPs = traffic_scaling * M_EPs


# print("uniform demand matrix:")

# M_EPs = gl.generate_uniform_traffic_demand_matrix(V, EPR)
# remote_link_flows, local_link_flows = _network.distribute_M_EPs_on_weighted_paths(ECMP_ASP, EPR, M_EPs)
# max_remote_link_load = np.max(remote_link_flows)/Cap_remote
# max_local_link_load = np.max(local_link_flows)/Cap_local
# # adapt the traffic scaling factor to 10x saturation
# traffic_scaling = 10.0/max(max_local_link_load, max_remote_link_load)
# M_EPs = traffic_scaling * M_EPs


print("starting nexullance_IT_fast test...")

nexu_it_fast = Nexullance_IT_fast_interface(V, arcs, 10.0, 10.0, False)
nexu_it_fast._initialize()
nexu_it_fast.set_parameters(V*2, 1000000, 0.00001)
nexu_it_fast_result = nexu_it_fast.add_next_matrix(M_EPs, 5, True)

print("====== nexu_IT_fast ========")
print("elasped_time [s]: ", nexu_it_fast_result.get_elapsed_time())
print("max_load: ", nexu_it_fast_result.get_max_core_link_load())
print("phi: ", nexu_it_fast_result.get_phi())
print("num_attempts: ", nexu_it_fast_result.get_num_attempts())

# # ==========================================================================
