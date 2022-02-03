#!/usr/bin/env python

import errno
import sys
from typing import Mapping, Tuple


def map_positions_to_alleles(locus_data_file: str) -> Mapping[int, Tuple[str, str]]:

    """
    Maps the position of a SNP to the alleles that were mapped.
    :param locus_data_file: A CSV file that contains the SNPs at the locus and their major/minor alleles.
    :return: A mapping from integer positions to a tuple of both minor and major allele.
    """

    CSV_POSITION_INDEX = 1
    CSV_MINOR_ALLELE_INDEX = 5
    CSV_MAJOR_ALLELE_INDEX = 6

    position_to_alleles = dict()

    with open(locus_data_file, 'r') as f_in:

        f_in.readline()  # Skip the header of the CSV

        # For each line, split by the delimiter ',' and extract relevant information
        for line in f_in:
            line = line.strip()
            if line:
                line_values = line.split(',')
                position_to_alleles[int(line_values[CSV_POSITION_INDEX])] = (
                    line_values[CSV_MINOR_ALLELE_INDEX], line_values[CSV_MAJOR_ALLELE_INDEX]
                )

    return position_to_alleles


def main():

    locus_data_file = sys.argv[1]
    reference_major_allele_mismatch_file = sys.argv[2]

    position_to_alleles = map_positions_to_alleles(locus_data_file)

    VCF_POSITION_INDEX = 1
    VCF_REF_ALLELE_INDEX = 3
    VCF_ALT_ALLELE_INDEX = 4

    from signal import signal, SIGPIPE, SIG_DFL
    signal(SIGPIPE, SIG_DFL)

    # Write out any SNPs where the reference allele does not match the major allele
    mismatches_out = open(reference_major_allele_mismatch_file, 'w')
    mismatches_out.write('SNP_Position\tReference_Allele\tAlternate_Allele\tMinor_Allele\tMajor_Allele\n')

    # Iterate over standard input stream
    for line in sys.stdin:

        line = line.rstrip()

        # Print all comments from the VCF File
        if line[0] == '#':
            print(line)
            continue

        # Split VCF file line by delimiter '\t'
        line_values = line.split('\t')

        # Extract the variant position, reference allele, and alternate allele from line
        position = int(line_values[VCF_POSITION_INDEX])
        ref_allele = line_values[VCF_REF_ALLELE_INDEX]
        alt_allele = line_values[VCF_ALT_ALLELE_INDEX]

        # Skip line if VCF position is not present in locus SNP data
        if position not in position_to_alleles:
            continue

        # Skip line if reference allele in VCF is not one of the locus alleles
        if ref_allele not in position_to_alleles[position]:
            continue

        # Skip line if alternate allele in VCF is not one of the locus alleles
        if alt_allele not in position_to_alleles[position]:
            continue
        
        # Report SNP if the reference allele does not match the major allele
        minor_allele, major_allele = position_to_alleles[position]
        if ref_allele != major_allele:
            mismatches_out.write(f'{position}\t{ref_allele}\t{alt_allele}\t{minor_allele}\t{major_allele}\n')
        
        # Output line to standard output
        print(line)

    mismatches_out.close()


if __name__ == '__main__':

    main()
