#!/usr/bin/env python

import argparse
from asyncore import read
import os
import re
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
    parser.add_argument('cell_type', help='Cell type of peak shapes to retrieve')
    parser.add_argument('output', help='Name of output file.')

    args = parser.parse_args()

    if not os.path.isfile(args.peaks):
        print(f'{args.peaks} is not a valid path to a file.', file=sys.stderr)
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


def read_peak_distributions(peaks, cell_type):

    '''
    Read peak distributions from read pile ups for the given cell type.

    :param peaks: A data frame containing peaks.
    :param cell_type: The cell type to retrieve read pile ups for.
    :return: A list of matrices of sample distributions for each peak.
    '''

    sample_dir = '/nfs/users/nfs_n/nm18/gains_team282/epigenetics/accessibility/merged/atac_seq/'

    samples = [sample_name for sample_name in os.listdir(sample_dir) if cell_type in sample_name]

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

    return samples, distributions


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


def main():

    args = parse_arguments()

    print('Reading Peaks from BED File')
    all_peaks = read_peaks(args.peaks)

    num_peak_sets = np.floor(len(all_peaks) / 200)
    peaks_sets = np.array_split(all_peaks, num_peak_sets)

    vectors_list = list()

    for i, peaks in enumerate(peaks_sets):

        print(f'Processing {i+1}/{num_peak_sets} Set of Peak')

        samples, distributions = read_peak_distributions(peaks, args.cell_type)

        vectors_list.append(fiedler_vectors(distributions))

    vectors = np.hstack(vectors_list)

    print('Saving to File')
    vectors_df = pd.DataFrame(vectors)
    vectors_df.columns = all_peaks.chr.str.cat(all_peaks.start.astype(str).str.cat(all_peaks.end.astype(str), sep='-'), sep=':').to_list()
    vectors_df.insert(0, 'Sample', samples)

    vectors_df.to_csv(args.output, index=False)


if __name__ == '__main__':

    main()
