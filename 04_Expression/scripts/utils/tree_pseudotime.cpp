#include <Rcpp.h>

#include <stdlib.h>

#include <cmath>
#include <iostream>
#include <vector>
#include <unordered_map>
#include <map>
#include <unordered_set>
#include <stack>

using namespace Rcpp;
using namespace std;

class DisjointSet {

private:

    unordered_map<size_t, size_t> parents;

public:

    explicit DisjointSet(const vector<size_t> &initial) {

        this->parents = unordered_map<size_t, size_t>();

        for (const size_t &value : initial) {
            this->parents[value] = value;
        }
    }

    size_t find_parent(const size_t &value) {

        if (this->parents[value] == value) {
            return value;
        }

        return this->find_parent(this->parents[value]);
    }

    void union_sets(const size_t &value_x, const size_t &value_y) {

        size_t parent_x = this->find_parent(value_x);
        size_t parent_y = this->find_parent(value_y);

        this->parents[parent_x] = parent_y;
    }
};

/**
 * Generates a minimum spanning tree (MST) of the projected points using Kruskal's algorithm. Given
 * a projected point that represents the root of the tree, calculates distances along the MST,
 * which are returned as the pseudotime for each sample.
 *
 * Let `k` be the number of centroids. Let `n` be the number of samples.
 *
 * @param distances A `n` by `n` matrix representing the distances between the projected points of
 * the samples.
 * @param root The index of sample in `projected` representing the root node.
 * @return A vector of `n` elements representing the pseudotime of the sample.
 */
// [[Rcpp::export]]
NumericVector tree_pseudotime(NumericMatrix distances, size_t root) {

    auto mst_adjacencies = unordered_map<size_t, vector<size_t>>();

    vector<size_t> disjoint_initial(distances.nrow());
    iota(disjoint_initial.begin(), disjoint_initial.end(), 0);
    DisjointSet disjoint(disjoint_initial);

    map<double, vector<size_t>> ordered_edges;
    for (size_t i = 0; i < distances.nrow(); i++) {
        for (size_t j = i; j < distances.nrow(); j++) {
            vector<size_t> edge = {i, j};
            ordered_edges[distances(i, j)] = edge;
        }
    }

    for (auto const &edge : ordered_edges) {

        size_t x = edge.second[0];
        size_t y = edge.second[1];

        if (disjoint.find_parent(x) != disjoint.find_parent(y)) {

            if (mst_adjacencies.count(x) == 0) {
                mst_adjacencies[x] = vector<size_t>();
            }

            if (mst_adjacencies.count(y) == 0) {
                mst_adjacencies[y] = vector<size_t>();
            }

            mst_adjacencies[x].push_back(y);
            mst_adjacencies[y].push_back(x);

            disjoint.union_sets(x, y);
        }
    }

    // Depth-First Traversal for Psuedotime

    auto pseudotime = NumericVector(distances.nrow());

    auto visited = unordered_set<size_t>();
    auto boundary = stack<size_t>();

    boundary.push(root);
    pseudotime[root] = 0;

    while (boundary.size() > 0) {

        size_t current = boundary.top();
        boundary.pop();
        visited.insert(current);

        for (const size_t &node : mst_adjacencies[current]) {
            if (visited.find(node) == visited.end()) {
                pseudotime[node] = distances(current, node) + pseudotime[current];
                boundary.push(node);
            }
        }
    }

    return pseudotime / max(pseudotime);
}
