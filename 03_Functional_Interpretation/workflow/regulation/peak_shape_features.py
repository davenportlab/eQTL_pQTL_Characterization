#!/usr/bin/env python

import argparse
import functools
import multiprocessing
import os
import sys

import pandas as pd
import numpy as np
import scipy.spatial.distance as scpd
import scipy.sparse.csgraph as sccs
import scipy.linalg as scl
import scipy.stats as scs
import pyBigWig


def parse_arguments():

    '''
    Parse and check command-line arguments.

    :return: Namespace containing command-line arguments.
    '''

    parser = argparse.ArgumentParser(description='Extract peak shape features for a cell type')

    parser.add_argument('peaks', help='Peaks in BED format')
    parser.add_argument('samples', help='File containing the list of samples to test')
    parser.add_argument('prefix', help='Prefix of output file.')
    parser.add_argument('--threads', help='Number of threads to use.', default='1')

    args = parser.parse_args()

    if not os.path.isfile(args.peaks):
        print(f'{args.peaks} is not a valid path to a file.', file=sys.stderr)
        exit(2)
    
    if not os.path.isfile(args.samples):
        print(f'{args.samples} is not a valid path to a file.', file=sys.stderr)
        exit(2)

    try:
        args.threads = int(args.threads)
    except ValueError:
        print(f'{args.threads} is not an integer.', file=sys.stderr)
        exit(2)
    
    return args


def read_peaks(peak_file):

    '''
    Read BED file containing peaks.

    :param peak_file: The BED file containing peaks.
    :return: A data frame containing peaks.
    '''

    peaks = pd.read_csv(
        peak_file, sep='\t',
        names=['chr', 'start', 'end'],
        usecols=[0, 1, 2],
        dtype={'chr': str, 'start': int, 'end': int}
    )

    return peaks


def read_peak_distributions(peaks, samples):

    '''
    Read peak distributions from read pile ups for the given cell type.

    :param peaks: A data frame containing peaks.
    :param samples_file: The samples to retrieve read pile ups for.
    :return: A list of matrices of sample distributions for each peak.
    '''

    sample_dir = '/nfs/users/nfs_n/nm18/gains_team282/epigenetics/accessibility/merged/atac_seq/'

    regions = [list() for _ in range(len(peaks))]

    for sample in samples:

        # For each sample, open BigWIG file for peak distributions
        bw_file = pyBigWig.open(os.path.join(sample_dir, sample, 'alignment', f'{sample}.bw'))

        for i, (_, row) in enumerate(peaks.iterrows()):

            # Retrieve BigWIG coverage in peak
            values = np.array(bw_file.values(row.chr, row.start, row.end))

            # Set any coverage that is not present (NaN) to zero
            values[np.isnan(values)] = 0

            # If the entire coverage is 0 (sample has no reads in this peak)
            # set to a uniform distribution. Otherwise, normalize coverage to
            # create a probability distribution.
            if values.sum() > 0:
                values = values / values.sum()
            else:
                values = np.array([1 / (row.end - row.start)] * (row.end - row.start))

            regions[i].append(values)

    distributions = [np.vstack(region) for region in regions]

    return distributions


def hellinger(x, y):

    '''
    Returns the Hellinger distance between two discrete probability distributions with the
    same support.

    :param x: A vector of probabilities representing a discrete probability distribution.
    :param y: A vector of probabilities representing a discrete probability distribution.
    :return: The Hellinger distance between the two distributions.
    '''

    return (1 / np.sqrt(2)) * np.sqrt(np.power(np.sqrt(x) - np.sqrt(y), 2).sum())


def fiedler_vectors(distributions):

    '''
    Retrieve the Fiedler vector for each peak.

    :param distributions: The list of peak distributions.
    :return: A matrix of Fiedler vectors.
    '''

    vectors = list()

    for distribution in distributions:

        distances_condensed = scpd.pdist(distribution, hellinger)

        distances = scpd.squareform(distances_condensed)

        laplacian = sccs.laplacian(1 - distances)

        eigenvalues, eigenvectors = scl.eig(laplacian)

        sort_indices = np.argsort(eigenvalues)

        eigenvectors = eigenvectors[:, sort_indices]

        vectors.append(np.real(eigenvectors[:, 1]))

    return np.transpose(np.vstack(vectors))


def coverage_pcs(distributions):

    '''
    Retrieve principal components for each distribution's coverage profile.

    :param distributions: The list of peak distributions.
    :return: A matrix of principal components.
    '''

    vectors = list()
    variance_explained = list()

    for distribution in distributions:

        coverage_stds = distribution.std(axis=0)
        distribution = distribution[:, np.flatnonzero(coverage_stds)]

        scaled = (distribution - distribution.mean(axis=0)) / distribution.std(axis=0)

        covariance = np.cov(scaled)

        eigenvalues, eigenvectors = np.linalg.eig(covariance)
        eigenvalues = np.real(eigenvalues)
        eigenvectors = np.real(eigenvectors)
        eigenvectors = eigenvectors[:, np.argsort(eigenvalues)[::-1]]

        vectors.append(eigenvectors[:, 0])

        variance_explained.append(eigenvalues[0] / np.sum(eigenvalues))

    return np.transpose(np.vstack(vectors)), np.array(variance_explained)


def process_peak_set(peaks_info, num_peak_sets, samples):

    i, peaks = peaks_info

    print(f'Processing Peak Set {i+1}/{num_peak_sets}')

    distributions = read_peak_distributions(peaks, samples)

    return fiedler_vectors(distributions), coverage_pcs(distributions)


def main():

    args = parse_arguments()

    print('Reading Peaks from BED File')
    all_peaks = read_peaks(args.peaks)

    num_peak_sets = np.floor(len(all_peaks) / 200)
    peaks_sets = np.array_split(all_peaks, num_peak_sets)

    print('Reading Samples')

    with open(args.samples, 'r') as f_in:

        samples = [sample_name.strip() for sample_name in f_in.readlines()]

    vectors_list = list()

    with multiprocessing.Pool(args.threads) as pool:

        vectors_list = pool.map(
            functools.partial(process_peak_set, num_peak_sets=num_peak_sets, samples=samples), 
            enumerate(peaks_sets)
        )

    fiedler_vector_matrix = np.hstack([items[0] for items in vectors_list])

    coverage_pc_matrix = np.hstack([items[1][0] for items in vectors_list])
    coverage_pc_variance_explained = np.concatenate([items[1][1] for items in vectors_list])

    print('Saving to File')

    fiedler_vector_df = pd.DataFrame(fiedler_vector_matrix)
    fiedler_vector_df.columns = all_peaks.chr.str.cat(all_peaks.start.astype(str).str.cat(all_peaks.end.astype(str), sep='-'), sep=':').to_list()
    fiedler_vector_df.insert(0, 'Sample', samples)

    fiedler_vector_df.to_csv(f'{args.prefix}_fiedler_vectors.csv', index=False)

    coverage_pc_df = pd.DataFrame(coverage_pc_matrix)
    coverage_pc_df.columns = all_peaks.chr.str.cat(all_peaks.start.astype(str).str.cat(all_peaks.end.astype(str), sep='-'), sep=':').to_list()
    coverage_pc_df.insert(0, 'Sample', samples)

    coverage_pc_df.to_csv(f'{args.prefix}_coverage_pcs.csv', index=False)

    coverage_pc_var_df = pd.DataFrame(coverage_pc_variance_explained)
    coverage_pc_var_df.columns = ['Variance_Explained']
    coverage_pc_var_df.insert(0, 'Peak', all_peaks.chr.str.cat(all_peaks.start.astype(str).str.cat(all_peaks.end.astype(str), sep='-'), sep=':'))

    coverage_pc_var_df.to_csv(f'{args.prefix}_coverage_pcs_var_explained.csv', index=False)


if __name__ == '__main__':

    main()
