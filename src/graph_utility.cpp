#include "graph_utility.hpp"
#include <boost/graph/dijkstra_shortest_paths_no_color_map.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/property_map/dynamic_property_map.hpp>  // For dynamic_properties in constructing a graph
#include <boost/graph/depth_first_search.hpp>

#include <queue>
#include <iostream>

Graph read_graph_from_arcs(int V, Eigen::MatrixX2i arcs, bool print) {
    // the input is a matrix (np array of size (E, 2)) that describes the arcs of the graph

    #ifdef DEBUG
    assert(arcs.rows() > 0 && "the matrix is empty");
    #endif

    // Create a graph object
    Graph g;
    // Setup dynamic properties to read vertex and edge properties
    dynamic_properties dp;
    dp.property("vertex_index", get(vertex_index, g));  // Read vertex indices

    // Create the vertices
    for (int i = 0; i < V; i++) {
        add_vertex(g);
    }

    // Create the edges
    for (int i = 0; i < arcs.rows(); i++) {
        int source = arcs(i, 0);
        int target = arcs(i, 1);
        add_edge(source, target, g);
    }    

    // Assign all edge weights as 1.0
    auto edge_pair = edges(g);
    for (auto it = edge_pair.first; it != edge_pair.second; ++it) {
        put(boost::edge_weight, g, *it, 1.0f);
    }

    if (print){
        std::cout << "Graph contents:" << std::endl;
        // Get the range of vertices
        auto vertexRange = boost::vertices(g);

        // Iterate over all vertices and print them
        std::cout << "Vertices:" << std::endl;
        for (auto it = vertexRange.first; it != vertexRange.second; ++it) {
            std::cout << *it << std::endl;
        }

        // Get the range of edges
        auto edgeRange = boost::edges(g);

        // Print out all edges in the format (int, int)
        std::cout << "Edges:" << std::endl;
        for (auto it = edgeRange.first; it != edgeRange.second; ++it) {
            auto source = boost::source(*it, g);
            auto target = boost::target(*it, g);
            std::cout << "(" << source << ", " << target << ")" << std::endl;
        }

        // Iterate through the edges and print their source, destination, and weight
        for (auto it = edge_pair.first; it != edge_pair.second; ++it) {
            auto src = source(*it, g);
            auto dest = target(*it, g);
            auto weight = get(boost::edge_weight, g, *it);
            std::cout << "Edge from vertex " << src << " to vertex " << dest
                    << " with weight " << weight << std::endl;
    }

    }
    return g;
}

// float compute_network_data_throughput(const float** M_R, float max_core_load, float access_link_cap){
//     // unit of flow: GBps

//     int n = sizeof(M_R) / sizeof(M_R[0]);
//     float sum_row = 0.0f;
//     float sum_col = 0.0f;
//     float sum_matrix = 0.0f;
//     float max_access_flow = 0.0f;

//     // first compute the sum of the rows, compare with max
//     for (int i = 0; i < n; i++) {
//         sum_row = 0.0f;
//         for (int j = 0; j < n; j++) {
//             sum_row += M_R[i][j];
//             if(sum_row > max_access_flow)
//                 max_access_flow = sum_row;
//         }
//         sum_matrix += sum_row;
//     }
//     // then compute the sum of the columns, compare with max
//     for (int i = 0; i < n; i++) {
//         sum_col = 0.0f;
//         for (int j = 0; j < n; j++) {
//             sum_col += M_R[j][i];
//             if(sum_col > max_access_flow)
//                 max_access_flow = sum_col;
//         }
//     }

//     float max_access_load = max_access_flow/access_link_cap;

//     if (max_access_load > max_core_load){
//         return sum_matrix/max_access_load;
//     }else{
//         return sum_matrix/max_core_load;
//     }

// }

std::pair<float, float> procress_M_EPs(const float** _M_EPs, const int num_vertices, const int _EPR, float*** out_M_R){


    int num_EPs = num_vertices*_EPR;

    // Calculate the max_access_link_load for later usage
    float maxEPflow = 0.0f; // maxRowColSum
    float total_flow = 0.0f;  // elementWiseSum

    // Array to store the sum of each column
    float* colSums = new float[num_EPs]{0};
    // Calculate the sums
    for (int i = 0; i < num_EPs; ++i) {
        float rowSum = 0.0f;
        for (int j = 0; j < num_EPs; ++j) {
            rowSum += _M_EPs[i][j];
            colSums[j] += _M_EPs[i][j];
            total_flow += _M_EPs[i][j];
        }
        // Update the maximum row sum
        maxEPflow = std::max(maxEPflow, rowSum);
    }
    // Check the column sums for the maximum
    for (int j = 0; j < num_EPs; ++j) {
        maxEPflow = std::max(maxEPflow, colSums[j]);
    }
    // Clean up the allocated memory
    delete[] colSums;

    float** temp_M_R = new float*[num_vertices];    // will be deleted in the caller
    for (int i = 0; i < num_vertices; i++) {
        temp_M_R[i] = new float[num_vertices]{0};    // will be deleted in the caller
    }
    for (int i = 0; i < num_vertices; i++) {
        for (int j = 0; j < num_vertices; j++) {
            if (i == j) {
                continue;
            }
            for (int sep = 0; sep < _EPR; sep++) {
                for (int dep = 0; dep < _EPR; dep++) {
                    temp_M_R[i][j] += _M_EPs[sep + i * _EPR][dep + j * _EPR];
                }
            }
        }
    }
    *out_M_R = temp_M_R;

    return std::make_pair(maxEPflow, total_flow);
}

void compute_all_shortest_paths_all_s_d(const Graph &G, std::vector<std::vector<Vertex>>** all_s_d_paths,
                                        const boost::property_map<Graph, boost::edge_weight_t>::type &weightmap){

        // iterate over all sources in the boost graph G
        for(size_t src = 0; src < num_vertices(G); src++){
            // dijkstra_shortest_paths_no_color_map(G, src, weight_map(weightmap));
            compute_all_shortest_paths_single_source(G, src,  all_s_d_paths[src], weightmap);
        }

}

void compute_all_shortest_paths_single_source(const Graph &G, Vertex s, std::vector<std::vector<Vertex>>* all_paths,
                                const boost::property_map<Graph, boost::edge_weight_t>::type &weightmap) {
    
    float dist[num_vertices(G)];
    boost::dijkstra_shortest_paths_no_color_map(G, s, distance_map(boost::make_iterator_property_map(dist, get(boost::vertex_index, G))).weight_map(weightmap));
    
    // Initialize paths from the source vertex
    all_paths[s].push_back({s});
    
    std::queue<Vertex> q;
    q.push(s);
    
    while (!q.empty()) {
        Vertex u = q.front();
        q.pop();

        OutEdgeIterator ei, ei_end;
        for (boost::tie(ei, ei_end) = boost::out_edges(u, G); ei != ei_end; ++ei) {
            Edge e = *ei;
            Vertex v = boost::target(e, G);

            if (dist[u] + weightmap[e] == dist[v]) {
                if (all_paths[v].empty()) {
                    q.push(v);
                }
                for (const auto &path : all_paths[u]) {
                    std::vector<Vertex> new_path = path;
                    new_path.push_back(v);
                    all_paths[v].push_back(new_path);
                }
            }
        }
    }
}


void compute_all_shortest_paths_single_s_d(const Graph &G, Vertex s, Vertex d, std::vector<std::vector<Vertex>> &result_all_paths,
                                const boost::property_map<Graph, boost::edge_weight_t>::type &weightmap) {
    std::vector<std::vector<Vertex>> all_paths[num_vertices(G)];
    
    #ifdef DEBUG
    assert(s != d);
    // make sure s and d are valid vertices in the graph
    assert(s < num_vertices(G) && d < num_vertices(G));
    #endif

    float dist[num_vertices(G)];
    boost::dijkstra_shortest_paths_no_color_map(G, s, distance_map(boost::make_iterator_property_map(dist, get(boost::vertex_index, G))).weight_map(weightmap));
    
    // Initialize paths from the source vertex
    all_paths[s].push_back({s});
    std::queue<Vertex> q;
    q.push(s);
    while (!q.empty()) {
        Vertex u = q.front();
        q.pop();

        OutEdgeIterator ei, ei_end;
        for (boost::tie(ei, ei_end) = boost::out_edges(u, G); ei != ei_end; ++ei) {
            Edge e = *ei;
            Vertex v = boost::target(e, G);

            if (dist[u] + weightmap[e] == dist[v]) {
                if (all_paths[v].empty()) {
                    q.push(v);
                }
                for (const auto &path : all_paths[u]) {
                    std::vector<Vertex> new_path = path;
                    new_path.push_back(v);
                    all_paths[v].push_back(new_path);
                }
            }
        }
        if (u == d){
            // properly copy the data from all_paths[u] to the result_all_paths
            result_all_paths = all_paths[u]; // This should be a deep copy
            return;
        }
    }
    // if we reach here, then there is no path from s to d
    // raise an error
    std::cout << "No path found from " << s << " to " << d << std::endl;
    #ifdef DEBUG
    assert(false);
    #endif
}

// // typedef typename property_traits<ColorMap>::value_type ColorValue;
// typedef color_traits<ColorValue> Color;


// // Define a custom DFS visitor for finding all paths from a source to all destinations with max length
// class all_destinations_path_finder_visitor : public boost::default_dfs_visitor {
// private:
//     Vertex source;
//     size_t max_length;
//     std::vector<Vertex>& current_path;
//     std::vector<std::vector<std::vector<std::vector<Vertex>>>>& all_paths_matrix;
//     std::vector<bool>& visited;
    
// public:
//     all_destinations_path_finder_visitor(Vertex s, size_t max_len, 
//                         std::vector<Vertex>& path, 
//                         std::vector<std::vector<std::vector<std::vector<Vertex>>>>& paths_matrix,
//                         std::vector<bool>& vis)
//         : source(s), max_length(max_len),
//           current_path(path), all_paths_matrix(paths_matrix), visited(vis) {
//         // Initialize path with source
//         current_path.clear();
//         current_path.push_back(source);
//     }
    
//     // Called when a vertex is discovered
//     template <typename Vertex, typename Graph>
//     void discover_vertex(Vertex u, const Graph& g) {        
//         // Add vertex to current path if it's not the source (already added)
//         if (u != source) {
//             current_path.push_back(u);
            
//             // If this is not the source and path is within max length, save the path
//             if (current_path.size() <= max_length + 1) { // +1 because path includes vertices, not just edges
//                 all_paths_matrix[source][u].push_back(current_path); // adding the current path

//             }
//         }
//     }
    
//     // Called when a vertex is finished (all neighbors processed)
//     template <typename Vertex, typename Graph>
//     void finish_vertex(Vertex u, const Graph& g) {
//         // TODO: color it white
//         boost::put(boost::vertex_color, g, u, boost::white_color);
//         visited[u] = false; // mark vertex as unvisited for future paths // TODO: color it white
        
//         // Backtrack: remove this vertex from the path
//         if (!current_path.empty() && current_path.back() == u) {
//             current_path.pop_back();
//         }
//     }
// };


// // this is a version that tries to call boost::depth_first_search with a custom visitor to find all paths from a source to all destinations with max length
// void compute_all_paths_with_max_length_all_s_d(const Graph &G, size_t max_length, 
//                                              std::vector<std::vector<std::vector<std::vector<Vertex>>>> &all_paths_matrix) {
//     size_t num_vertices = boost::num_vertices(G);
    
//     // Initialize the all_paths_matrix with the right dimensions
//     all_paths_matrix.resize(num_vertices);
//     for (size_t src = 0; src < num_vertices; ++src) {
//         all_paths_matrix[src].resize(num_vertices);
//     }
    
//     // Find paths for all sources to all destinations
//     for (size_t src = 0; src < num_vertices; ++src) {

//         std::vector<Vertex> current_path;
//         std::vector<bool> visited(num_vertices, false);
//         // Create a color map (required by boost DFS)
//         std::vector<boost::default_color_type> color_map(num_vertices, boost::white_color);
//         // Create the visitor to find paths from src to all destinations
//         all_destinations_path_finder_visitor vis(src, max_length, current_path, all_paths_matrix, visited);
//         // Run a custom DFS using our visitor
//         depth_first_search(G, 
//                            boost::visitor(vis)
//                            .edge_filter(boost::ref(vis))
//                            .color_map(boost::make_iterator_property_map(color_map.begin(), get(boost::vertex_index, G)))
//                            .root_vertex(src)
//                           );
        
//         // Check if we found paths to all destinations
//         for (size_t dst = 0; dst < num_vertices; ++dst) {
//             if (src == dst) continue; // Skip self-loops
            
//             // If no paths were found, log a message
//             if (all_paths_matrix[src][dst].empty()) {
//                 std::cout << "No paths found from " << src << " to " << dst 
//                           << " with length <= " << max_length << std::endl;
//                 #ifdef DEBUG
//                 assert(false);
//                 #endif
//             }
//         }
//     }
// }


typedef std::pair<Vertex, std::pair<boost::optional<Edge>, std::pair<OutEdgeIterator, OutEdgeIterator> > > VertexInfo;

// a variation of depth first search, simple but maybe slow
// If no paths were found between any src and dest vertices, log a message
void compute_all_paths_with_max_length_all_s_d(const Graph &g, size_t max_length, 
                                             std::vector<std::vector<std::vector<std::vector<Vertex>>>> &result_all_paths) {
    size_t num_vertices = boost::num_vertices(g);

    // Initialize the result_all_paths with the right dimensions
    result_all_paths.resize(num_vertices);
    for (size_t src = 0; src < num_vertices; ++src) {
        result_all_paths[src].resize(num_vertices);
    }

    // iterate over all sources in the boost graph g
    VertexIterator v_iter, v_iter_end;
    Vertex u;
    for (boost::tie(v_iter, v_iter_end) = boost::vertices(g); v_iter != v_iter_end; ++v_iter) {

        std::vector<Vertex> current_path; // For maintaining insertion order
        current_path.reserve(max_length + 1); // Reserve space for max length + 1 vertices
        std::unordered_set<Vertex> current_path_set; // For O(1) lookups

        boost::optional<Edge> src_e;
        OutEdgeIterator ei, ei_end;
        std::vector<VertexInfo> stack;

        // Possible optimization for vector
        //stack.reserve(num_vertices(g));

        current_path.clear();
        current_path_set.clear();
        current_path.push_back(*v_iter);
        current_path_set.insert(*v_iter);
        boost::tie(ei, ei_end) = out_edges(*v_iter, g);
        stack.push_back(std::make_pair(*v_iter, std::make_pair(boost::optional<Edge>(), std::make_pair(ei, ei_end))));

        while (!stack.empty()) {
        VertexInfo& back = stack.back();
        u = back.first;
        src_e = back.second.first;
        boost::tie(ei, ei_end) = back.second.second;
        stack.pop_back();
        while (ei != ei_end) {
            Vertex v = target(*ei, g);
            // Check if path length exceeds max_length
            if ((current_path.size() < max_length + 1) &&                    // if the current path length is ok
                (current_path_set.find(v) == current_path_set.end())) {      // if v is not already in the current path - O(1) lookup
                // vis.tree_edge(*ei, g);
                src_e = *ei;
                stack.push_back(std::make_pair(u, std::make_pair(src_e, std::make_pair(++ei, ei_end))));
                u = v;
                current_path.push_back(u);
                current_path_set.insert(u);

                // Add the current path to the result_all_paths, 
                // assert that the path does not already exists
                #ifdef DEBUG
                for (const auto& path : result_all_paths[*v_iter][u]) {
                    assert(path != current_path && "Path already exists in result_all_paths");
                }
                #endif

                result_all_paths[*v_iter][u].push_back(current_path); // adding the current path
                #ifdef DEBUG
                assert(current_path.size() <= max_length + 1 && "Path exceeds max length");
                std::cout << "Found path from " << *v_iter << " to " << u << ": ";
                for (const auto& vertex : current_path) std::cout << vertex << " " ;
                std::cout << std::endl;
                #endif
                boost::tie(ei, ei_end) = out_edges(u, g);
            } else {
                ++ei;
            }
        }
        #ifdef DEBUG
        assert(!current_path.empty() && current_path.back() == u);
        #endif
        current_path_set.erase(u); // Remove from set before popping from vector
        current_path.pop_back(); // Backtrack: remove this vertex from the path
        }      
    }

    #ifdef DEBUG
    // assert that there exists at least one path from each source to each destination
    for (size_t src = 0; src < num_vertices; ++src) {
        for (size_t dst = 0; dst < num_vertices; ++dst) {
            if (src == dst) continue; // Skip self-loops
            if (result_all_paths[src][dst].empty()) {
                std::cout << "Exception: No paths found from " << src << " to " << dst
                          << " with length <= " << max_length << std::endl;
                assert(false);
            }
        }
    }
    #endif
}
