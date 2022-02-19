#include <Rcpp.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <vector>

using namespace Rcpp;
using namespace std;

/**
 * Calculates the Context Likelihood of Relatedness (CLR) for all gene pairs in an adjacency
 * matrix.
 *
 * @param adjacency A square adjacency matrix.
 * @return A transformed adjacency matrix with CLR values.
 */
// [[Rcpp::export]]
NumericMatrix clr(NumericMatrix adjacency) {

    std::vector<double> adjacency_totals;

    for (size_t i = 0; i < adjacency.nrow(); i++) {

        double total = 0;
        for (size_t j = 0; j < adjacency.ncol(); j++) {

            total += adjacency(i, j);
        }

        adjacency_totals.push_back(total);
    }

    std::vector<double> adjacency_means;

    for (size_t i = 0; i < adjacency.nrow(); i++) {

        double mean = adjacency_totals[i] / adjacency.ncol();
        adjacency_means.push_back(mean);
    }

    std::vector<double> adjacency_sd;

    for (size_t i = 0; i < adjacency.nrow(); i++) {

        double deviations = 0;
        for (size_t j = 0; j < adjacency.ncol(); j++) {

            deviations += pow(adjacency(i, j) - adjacency_means[i], 2);
        }

        double sd = sqrt(deviations / (adjacency.ncol() - 1));
        adjacency_sd.push_back(sd);
    }
    
    // Z[i,j] is the z-score of gene j in the marginal distribution of gene i
    NumericMatrix z_scores(adjacency.nrow(), adjacency.ncol());

    for (size_t i = 0; i < adjacency.nrow(); i++) {

        for (size_t j = 0; j < adjacency.ncol(); j++) {

            z_scores(i, j) = (adjacency(i, j) - adjacency_means[i]) / adjacency_sd[i];
        }
    }

    NumericMatrix clr(z_scores.nrow(), z_scores.ncol());

    for (size_t i = 0; i < z_scores.nrow(); i++) {

        for (size_t j = 0; j < z_scores.ncol(); j++) {

            double z_i = z_scores(j, i);
            double z_j = z_scores(i, j);
            clr(i, j) = sqrt((z_i * z_i) + (z_j * z_j));
        }
    }

    return clr;
}
