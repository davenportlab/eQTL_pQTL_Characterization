#----------------------------------------------------------
# Fine Mapping with FINEMAP
# Created: 18 February 2022
#----------------------------------------------------------

#----------------------------------------------------------
# Environment
#----------------------------------------------------------

import os
import re
import sys

import numpy as np
import pandas as pd

cis_or_trans = sys.argv[1]

#----------------------------------------------------------
# Iterate over Credible Sets
#----------------------------------------------------------

cred_output = {
    'Gene': list(),
    'k': list(),
    'Post_Prob_k': list(),
    'Credibility_Set': list(),
    'Credibility_Set_Min_LD': list(),
    'Credibility_Set_Mean_LD': list(),
    'Credibility_Set_Median_LD': list(),
    'SNP': list(),
    'Prob_SNP_in_CS': list()
}

for file_name in os.listdir('credible_sets/'):

    file_path = f'credible_sets/{file_name}'

    gene_name = re.sub('\\.cred.*', '', file_name)

    cred_min_ld = list()
    cred_mean_ld = list()
    cred_median_ld = list()

    with open(file_path, 'r') as file_in:
        
        first_line = file_in.readline().strip()
        k = int(re.sub('\).*', '', re.sub('.* is ', '', first_line)))
        k_post_prob = float(re.sub('.*= ', '', first_line))

        file_in.readline()  # Skip second line containing log(Bayes Factor)

        min_ld_line = file_in.readline().strip()
        mean_ld_line = file_in.readline().strip()
        median_ld_line = file_in.readline().strip()

        cred_min_ld = [float(x) for i, x in enumerate(re.sub('.*\) ', '', min_ld_line).split()) if i % 2 == 0]
        cred_mean_ld = [float(x) for i, x in enumerate(re.sub('.*\) ', '', mean_ld_line).split()) if i % 2 == 0]
        cred_median_ld = [float(x) for i, x in enumerate(re.sub('.*\) ', '', median_ld_line).split()) if i % 2 == 0]

    cred = pd.read_csv(file_path, sep=' ', skiprows=5, index_col=0)

    for signal in range(k):

        signal_cred = cred.iloc[:, (2 * signal):(2 * signal + 2)].dropna()
        n_snps = signal_cred.shape[0]

        cred_output['Gene'] += [gene_name] * n_snps
        cred_output['k'] += [k] * n_snps
        cred_output['Post_Prob_k'] += [k_post_prob] * n_snps
        cred_output['Credibility_Set'] += [signal + 1] * n_snps
        cred_output['Credibility_Set_Min_LD'] += [cred_min_ld[signal]] * n_snps
        cred_output['Credibility_Set_Mean_LD'] += [cred_mean_ld[signal]] * n_snps
        cred_output['Credibility_Set_Median_LD'] += [cred_median_ld[signal]] * n_snps
        cred_output['SNP'] += signal_cred.iloc[:, 0].tolist()
        cred_output['Prob_SNP_in_CS'] += signal_cred.iloc[:, 1].tolist()


cred_output_df = pd.DataFrame(cred_output).sort_values(
    by=['Gene', 'Post_Prob_k', 'Credibility_Set', 'Prob_SNP_in_CS'], 
    ascending=[True, False, True, False]
)

cred_output_df.to_csv(f'{cis_or_trans}_credible_sets.tsv', sep='\t', index=False)

#----------------------------------------------------------
# Iterate over SNP Data
#----------------------------------------------------------

# Create ouptut file with header from one of the input files
file_name = os.listdir('snps/')[0]
file_path = f'snps/{file_name}'
snps = pd.read_csv(file_path, sep=' ', index_col=0)

with open(f'{cis_or_trans}_pips.tsv', 'w') as f_out:
    header = snps.columns.tolist()
    header = ['gene'] + header + ['\n']
    f_out.write('\t'.join(header))

for file_name in sorted(os.listdir('snps/')):

    file_path = f'snps/{file_name}'

    gene_name = re.sub('\\.snp', '', file_name)

    snps = pd.read_csv(file_path, sep=' ', index_col=0).sort_values(by='position')
    snps.insert(0, 'Gene', [gene_name] * snps.shape[0])

    snps.to_csv(f'{cis_or_trans}_pips.tsv', sep='\t', header=False, index=False, mode='a')
