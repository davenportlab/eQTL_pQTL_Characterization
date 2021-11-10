#include <Rcpp.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>

using namespace Rcpp;
using namespace std;

/**
 * Calculates the frequency of items in a vector. The vector must contain `n` items, with each
 * element in the vector representing the cluster label. The cluster labels must be from the set of
 * integers [0, n-1].
 * 
 * @param vec A vector with cluster assignments for `n` items into k=`clusters` partitions.
 * @param out_array A pointer to an array of doubles into which the frequency array is deposited.
 *    The array size is assumed to be `clusters`.
 * @param clusters The number of partitions of the `n` items.
 * @param n The number of items that were clustered.
 */
void frequency_array(NumericVector vec, double *out_array, int clusters, int n) {
  
  // Initialize array to 0
  for (int i = 0; i < clusters; i++) {
    out_array[i] = 0.0;
  }
  
  // Iterate over vector and update counts
  for (int i = 0; i < n; i++) {
    out_array[(int)vec[i]] += 1.0;
  }
  
  // Divide by total count
  for (int i = 0; i < clusters; i++) {
    out_array[i] /= (double)n;
  }
}

/**
 * Calculates the frequency matrix of the intersections of sets within two partitions. The vectors
 * `x` and `y` must contain `n` items each. Each element within the vectors must represent a
 * cluster label. The cluster labels must be from the set of integers [0, n-1].
 * 
 * @param x A vector with cluster assignments for `n` items into k=`n_x` partitions.
 * @param y A vector with cluster assignments for `n` items into k=`n_y` partitions.
 * @param out_array A pointer to an array of arrays of doubles into which the frequency matrix is 
 *    deposited.
 * @param n_x The number of sets within the partition of the `n` items in the first clustering.
 * @param n_y The number of sets within the partition of the `n` items in the second clustering.
 * @param n The number of items that were clustered.
 */
void frequency_matrix(NumericVector x, NumericVector y, double **out_array, int n_x, int n_y, int n) {
  
  // Initialize matrix to 0
  for (int i = 0; i < n_x; i++) {
    for (int j = 0; j < n_y; j++) {
      out_array[i][j] = 0.0;
    }
  }
  
  // Iterate over vector and update intersection counts
  for (int i = 0; i < n; i++) {
    out_array[(int)x[i]][(int)y[i]] += 1.0;
  }
  
  // Divide by total count
  for (int i = 0; i < n_x; i++) {
    for (int j = 0; j < n_y; j++) {
      out_array[i][j] /= (double)n;
    }
  }
}

/**
 * Calculates the Variation of Information (VI) between two partitions of `n` items. The vectors
 * `x` and `y` must contain `n` items each. Each element within the vectors must represent a
 * cluster label. The cluster labels must be from the set of integers [0, n-1].
 * 
 * @param x A vector with cluster assignments for `n` items into k=`n_x` partitions.
 * @param y A vector with cluster assignments for `n` items into k=`n_y` partitions.
 * @param n_x The number of sets within the partition of the `n` items in the first clustering.
 * @param n_y The number of sets within the partition of the `n` items in the second clustering.
 */
// [[Rcpp::export]]
double vi(NumericVector x, int n_x, NumericVector y, int n_y) {
  
  // Assumes that x and y have the same length
  size_t n = x.size();
  
  // Initialize p and q vectors
  
  double *p = new double[n_x];
  double *q = new double[n_y];
  
  frequency_array(x, p, n_x, n);
  frequency_array(y, q, n_y, n);
  
  // Initialize R matrix
  
  double **r = (double**) malloc(n_x * sizeof(double*));
  for (int i = 0; i < n_x; i++) {
    r[i] = (double*) malloc(n_y * sizeof(double));
  }
  
  frequency_matrix(x, y, r, n_x, n_y, n);
  
  // Calculate Variation of Information
  
  double vi_result = 0.0;
  
  for (int i = 0; i < n_x; i++) {
    for (int j = 0; j < n_y; j++) {
      
      double r_ij = r[i][j];
      double p_i = p[i];
      double q_j = q[j];
      
      double log_p = r_ij > 0 && p_i > 0 ? log2(r_ij / p[i]) : 0.0;
      double log_q = r_ij > 0 && q_j > 0 ? log2(r_ij / q[j]) : 0.0;
      
      vi_result -= r_ij * (log_p + log_q);
    }
  }
  
  // Free Memory
  
  delete[] p;
  delete[] q;
  
  for (int i = 0; i < n_x; i++) {
    free(r[i]);
  }
  free(r);
  
  return vi_result;
}
