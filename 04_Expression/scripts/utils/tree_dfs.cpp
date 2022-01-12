#include <Rcpp.h>

#include <stdlib.h>

#include <math.h>

#include <iostream>
#include <unordered_map>
#include <vector>
#include <stack>
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
 * Performs Depth-First Search (DFS) on the `tree` to identify tree branches. Each centroid will 
 * be assigned a branch number from 0 up to `k-1`, depending on the number of branches detected.
 *
 * Let `k` be the number of centroids.
 *
 * @param tree A `k-1` by 2 matrix representing the adjacency list of the tree. The adjacency list
 * must use sample indices starting from 0. The matrix must be sorted by the first column.
 * @return A vector with `k` elements representing the branch assignments. Branches are assigned
 * from 0 up to `k-1`. A branch node is represented as a negative number, where the absolute value
 * represents the degree of the branch point.
 */
// [[Rcpp::export]]
NumericVector tree_dfs(NumericMatrix tree) {

    adj_map map;
    adjacency_matrix_to_map(tree, map);

    auto branches = NumericVector(tree.nrow() + 1);

    // Find the first node that has a degree of 1.
    //  We begin the DFS traversal from this leaf.

    int first_node = -1;
    for (auto &element : map) {
        if (element.second.size() == 1) {
            first_node = element.first;
            break;
        }
    }

    // Represents the current branch that the traversal is on.
    //  Each time the traversal hits a branch point, this is updated.
    int current_branch = 0;

    // The stack is used to identify the next node to traverse.
    // The set of visited nodes avoids backtracking over the graph.
    auto nodes = stack<int>();
    auto visited_nodes = unordered_set<int>();

    // Push the first leaf node onto the stack.
    nodes.push(first_node);

    // Set the branch of the first node.
    branches[first_node] = current_branch;

    // Iterate until the stack is empty (all nodes have been visited)
    while (nodes.size() > 0) {

        // Retrieve the node on top of the stack
        int node = nodes.top();
        nodes.pop();

        // Note down that the node has been visited
        visited_nodes.insert(node);

        // If the number of adjacent nodes is 2 or less, we are dealing with a node that is on a
        //  branch. In this case, we will continue with the same branch number.
        if (map[node].size() <= 2) {
            
            // For each adjacent node that has not yet been visited, add the node to the stack and
            //  set the branch of the adjacent node to the current branch.
            for (int adjacent_node : map[node]) {
                if (visited_nodes.find(adjacent_node) == visited_nodes.end()) {
                    
                    nodes.push(adjacent_node);
                    branches[adjacent_node] = branches[node];
                }
            }

        } 
        // Otherwise, we are dealing with a branch node.
        else {

            // If the current node is a branch node, set the branch to -1 * degree(node).
            branches[node] = -1 * (int)map[node].size();

            // For each adjacent node that has not yet been visited, add the node to the stack and
            //  set the branch of the adjacent node to a new branch.
            for (int adjacent_node : map[node]) {
                if (visited_nodes.find(adjacent_node) == visited_nodes.end()) {
                    
                    current_branch += 1;

                    nodes.push(adjacent_node);
                    branches[adjacent_node] = current_branch;
                }
            }
        }
    }

    return branches;
}
