#include <Rcpp.h>

#include <stdlib.h>

#include <math.h>

#include <iostream>
#include <unordered_map>
#include <vector>
#include <queue>
#include <unordered_set>

using namespace Rcpp;
using namespace std;

typedef unordered_map<int, vector<int>> adj_map;

/**
 * Converts the adjacency matrix of the tree to a dictionary for efficient access to endpoints.
 *
 * @param tree A `D` by 2 matrix representing the adjacency list of the tree.
 * @param map A mapping from the sample identifier (`int`) to a `vector<int>` of connected samples.
 */
void adjacency_matrix_to_map(NumericMatrix &tree, adj_map &map) {
    
    for (int i = 0; i < tree.nrow(); i++) {

        int node_1 = tree(i, 0);
        int node_2 = tree(i, 1);

        if (map.count(node_1) == 0) {
            map[node_1] = vector<int>();  
        }

        map[node_1].push_back(node_2);

        if (map.count(node_2) == 0) {
            map[node_2] = vector<int>();
        }

        map[node_2].push_back(node_1);
    }
}

/**
 * Performs a Breadth-First Search (BFS) on the adjacency map from the `start_node` to the
 * `end_node`. Outputs the nodes traversed into `path`.
 * 
 * @param tree The adjacency map representing the spanning tree.
 * @param start_node The node at which to begin the search.
 * @param end_node The node to search for within the tree.
 * @param latent The centroids in the latent space to traverse along.
 * @param traversed_path The vector to fill with nodes traversed on the path to the `end_node`.
 * @return The total distance required to traverse the path.
 */
double bfs(adj_map &tree, int start_node, int end_node, NumericMatrix &latent, vector<int> &traversed_path) {

    queue<vector<int>> paths;
    unordered_set<int> visited_nodes = unordered_set<int>();

    paths.push(vector<int> { start_node });
    visited_nodes.insert(start_node);

    while (paths.size() > 0) {

        vector<int> path = paths.front();
        paths.pop();

        int path_last_node = path.back();

        if (path_last_node == end_node) {
            traversed_path = path;
        }

        visited_nodes.insert(path_last_node);
        
        for (int adjacent_node : tree[path_last_node]) {
            if (visited_nodes.find(adjacent_node) == visited_nodes.end()) {
                vector<int> new_path = path;
                new_path.push_back(adjacent_node);
                paths.push(new_path);
            }
        }
    }

    double distance = 0.0;

    for (size_t i = 0; i < traversed_path.size() - 1; i++) {
        
        int node_i = traversed_path[i];
        int node_j = traversed_path[i + 1];

        // R uses 1-indexing, but C++ uses 0-indexing
        node_i -= 1;
        node_j -= 1;

        // Calculate distance between adjacent pairs of latent variables
        NumericVector node_i_latent_coord = latent(_, node_i);
        NumericVector node_j_latent_coord = latent(_, node_j);
        NumericVector between = node_j_latent_coord - node_i_latent_coord;

        distance += sqrt(sum(between * between));
    }

    return distance;
}

/**
 * Performs Breadth-First Search (BFS) on the `tree` to identify the path between samples of
 * interest in the latent space. For each pair of samples, the output will be a vector in the
 * latent space, tangent to the spanning tree, and pointing from the first sample to the second 
 * sample.
 *
 * Let `n` be the number of sample pairs being queried. Let `D` be the dimension of the input
 * space. Let `d` be the dimension of the latent space.
 *
 * @param tree A `D` by 2 matrix representing the adjacency list of the tree.
 * @param samples An `n` by 2 matrix representing pairs of samples.
 * @param latent A `d` by `D` matrix representing the sample centroids in the latent space.
 * @return A `n` by `d` matrix of `d`-dimensional vectors. The direction of the vector is aligned
 * with the path on the manifold between the two samples.
 */
// [[Rcpp::export]]
NumericMatrix stree_bfs(NumericMatrix tree, NumericMatrix samples, NumericMatrix latent) {

    adj_map map;
    adjacency_matrix_to_map(tree, map);

    NumericMatrix output = NumericMatrix(latent.nrow() + 1, samples.nrow());

    for (size_t i = 0; i < samples.nrow(); i++) {

        int start_node = samples(i, 0);
        int end_node = samples(i, 1);

        vector<int> path = {};

        double distance = bfs(map, start_node, end_node, latent, path);

        // Use the vector between the first and second nodes of the path
        //  as an approximation of the tangent to the path on the graph.
        int second_node = path[1];

        NumericVector start_node_coord = latent(_, start_node - 1);
        NumericVector second_node_coord = latent(_, second_node - 1);
        NumericVector between = second_node_coord - start_node_coord;

        between = between / sum(sqrt(between * between));

        size_t j;
        for (j = 0; j < latent.nrow(); j++) {
            output(j, i) = between[j];
        }
        output(j, i) = distance;
    }

    return output;
}
