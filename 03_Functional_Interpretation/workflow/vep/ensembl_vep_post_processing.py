#!/usr/bin/env python

import os
import subprocess


os.chdir('/nfs/users/nfs_n/nm18/gains_team282/epigenetics/vep/ensembl_vep/')

# Put header CSQ info into output file
header_terms = list()
with open('csq_format.txt') as header_file:
    header_terms = header_file.readline().strip().split('|')

vcf_snp_data = subprocess.run(['bcftools', 'query', '-f', '%CHROM\t%POS\t%ID\t%REF\t%ALT\n', 'qtl_both_ref.vep.vcf.gz'], capture_output=True, text=True)
vcf_snp_data = vcf_snp_data.stdout.split('\n')

vep_data = subprocess.run(['bcftools', 'query', '-f', '%INFO/CSQ\n', 'qtl_both_ref.vep.vcf.gz'], capture_output=True, text=True)
vep_data = [[x.replaceAll('|', '\t') for x in line.strip().split(',')] for line in vep_data.stdout.split('\n')]

with open('qtl_vep_output.tsv', 'w') as out_file:

    vep_terms = '\t'.join(header_terms)

    out_file.write(f'CHROM\tPOS\tID\tREF\tALT\t{vep_terms}\n')

    for snp_info, vep_info in zip(vcf_snp_data, vep_data):

        for consequence in vep_info:

            out_file.write(f'{snp_info}\t{consequence}\n')

