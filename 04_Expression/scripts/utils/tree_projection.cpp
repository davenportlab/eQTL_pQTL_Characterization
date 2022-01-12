#include <Rcpp.h>

#include <stdlib.h>

#include <cmath>

#include <iostream>

using namespace Rcpp;
using namespace std;

/**
 * Performs a projection of the latent points onto the tree. 
 *
 * Let `k` be the number of centroids. Let `n` be the number of samples.
 *
 * @param latent A 2 by `n` matrix representing the latent points of the samples.
 * @param centroids A 2 by `k` matrix representing the connected centroids of the graph.
 * @param tree A `k-1` by 2 matrix representing the adjacency list of the tree. The adjacency list
 * must use sample indices starting from 0. The matrix must be sorted by the first column.
 * @return A `n` by `4` matrix representing the projection of the latent points onto the tree. The
 * first two values are the coordinates of the point. The last two values are the centroids that
 * form the line onto which the latent point was projected.
 */
// [[Rcpp::export]]
NumericMatrix tree_projection(NumericMatrix latent, NumericMatrix centroids, NumericMatrix tree) {

    NumericMatrix projection(latent.ncol(), 4);

    // Iterate over all the latent points
    for (size_t i = 0; i < latent.ncol(); i++) {

        // For each point, keep track of the segment and projected point that minimize distance
        double min_dist = -1.0;
        size_t min_sample_i = 0;
        size_t min_sample_j = 0;
        double min_x = 0;
        double min_y = 0;

        // Iterate over all segments of the tree to identify the closest segment
        for (size_t j = 0; j < tree.nrow(); j++) {

            // Create references or extract information from the input matrices
            double &zx = latent(0, i);
            double &zy = latent(1, i);

            size_t sample_i = (size_t) tree(j, 0);
            size_t sample_j = (size_t) tree(j, 1);

            double &ax = centroids(0, sample_i);
            double &ay = centroids(1, sample_i);

            double &bx = centroids(0, sample_j);
            double &by = centroids(1, sample_j);

            // Calculate projection factor (indicates the scaling of (b - a) to reach the projected point from a
            double t = (((zx - ax) * (bx - ax)) + ((zy - by) * (by - ay))) / (((bx - ax) * (bx - ax)) + ((by - ay) * (by - ay)));

            // Calculate the distance from the latent point to the segment
            double dist = 0.0;
            if (0.0 <= t && t <= 1.0) {
                dist = abs(((zx - ax) * (zy - by)) - ((zx - bx) * (zy - ay))) / sqrt(((bx - ax) * (bx - ax)) + ((by - ay) * (by - ay)));
            }
            else {
                double dist_a = sqrt(((ax - zx) * (ax - zx)) + ((ay - zy) * (ay - zy)));
                double dist_b = sqrt(((bx - zx) * (bx - zx)) + ((by - zy) * (by - zy)));
                dist = fmin(dist_a, dist_b);
            }

            // Update if the distance is the new minimum distance
            if (min_dist < 0.0 || min_dist > dist) {
                min_dist = dist;
                min_sample_i = sample_i;
                min_sample_j = sample_j;
                min_x = ax + (bx - ax) * t;
                min_y = ay + (by - ay) * t;
            }
        }

        // Update projection matrix
        projection(i, 0) = min_x;
        projection(i, 1) = min_y;
        projection(i, 2) = min_sample_i;
        projection(i, 3) = min_sample_j;
    }

    return projection;
}
