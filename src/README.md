# Nexullance_IT_boost

## to make

first set correct python path and pybind11 in CMakeLists.txt
install eigen3, e.g., sudo apt install libeigen3-dev
install boost, e.g., sudo apt install libboost-all-dev

#### definitions.hpp: 
Defines boost graph data structures

#### graph_utility.hpp, graph_utility.cpp: 
Utility functions for graph algorithms, includes graph-construction (from a list of arcs) and path-finding algorithms.
Sep2024: Added a function to calculate the network data throughput (given demand matrix and max core link load)

#### Nexullance_IT.hpp, Nexullance_IT.cpp: 
Implementation of the Nexullance_IT algorithm

#### MD_Nexullance_IT.hpp , MD_Nexullance_IT.cpp: 
Implementation of the MD_Nexullance_IT algorithm
  
#### lib_nexullance_IT.hpp, lib_nexullance_IT.cpp: 
pybind11 interface, wraps around Nexullance_IT and MD_Nexullance_IT classes
 





