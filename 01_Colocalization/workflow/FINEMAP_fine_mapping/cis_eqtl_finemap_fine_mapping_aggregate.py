#----------------------------------------------------------
# Fine Mapping with FINEMAP
# Created: 18 February 2022
#----------------------------------------------------------

#----------------------------------------------------------
# Environment
#----------------------------------------------------------

import glob
import os
import re
import sys

import numpy as np
import pandas as pd

chr_input = sys.argv[1]

#----------------------------------------------------------
# Iterate over Credible Sets (Full)
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

for file_path in glob.iglob('full_cred_sets/*.cred*'):

    gene_name = re.sub('\\.cred.*', '', os.path.basename(file_path))

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

cred_output_df.to_csv(f'full_chr{chr_input}_credible_sets.tsv', sep='\t', index=False)

#----------------------------------------------------------
# Iterate over Credible Sets (Conditional)
#----------------------------------------------------------

cred_output = {
    'Gene': list(),
    'Signal': list(),
    'Credibility_Set_Min_LD': list(),
    'Credibility_Set_Mean_LD': list(),
    'Credibility_Set_Median_LD': list(),
    'SNP': list(),
    'Prob_SNP_in_CS': list()
}

for file_path in glob.iglob('conditional_cred_sets/*.cred*'):

    locus_name = re.sub('\\.cred.*', '', os.path.basename(file_path))
    gene_name = re.sub('-.*', '', locus_name)
    signal = re.sub('.*-', '', locus_name)

    cred_min_ld = list()
    cred_mean_ld = list()
    cred_median_ld = list()

    with open(file_path, 'r') as file_in:
        
        file_in.readline()  # Skip first line containing posterior probability for k
        file_in.readline()  # Skip second line containing log(Bayes Factor)

        min_ld_line = file_in.readline().strip()
        mean_ld_line = file_in.readline().strip()
        median_ld_line = file_in.readline().strip()

        cred_min_ld = float(re.sub('.*\) ', '', min_ld_line).split()[0])
        cred_mean_ld = float(re.sub('.*\) ', '', mean_ld_line).split()[0])
        cred_median_ld = float(re.sub('.*\) ', '', median_ld_line).split()[0])

    cred = pd.read_csv(file_path, sep=' ', skiprows=5, index_col=0)

    signal_cred = cred.iloc[:, 0:2].dropna()
    n_snps = signal_cred.shape[0]

    cred_output['Gene'] += [gene_name] * n_snps
    cred_output['Signal'] += [signal] * n_snps
    cred_output['Credibility_Set_Min_LD'] += [cred_min_ld] * n_snps
    cred_output['Credibility_Set_Mean_LD'] += [cred_mean_ld] * n_snps
    cred_output['Credibility_Set_Median_LD'] += [cred_median_ld] * n_snps
    cred_output['SNP'] += signal_cred.iloc[:, 0].tolist()
    cred_output['Prob_SNP_in_CS'] += signal_cred.iloc[:, 1].tolist()


cred_output_df = pd.DataFrame(cred_output).sort_values(
    by=['Gene', 'Signal', 'Prob_SNP_in_CS'], 
    ascending=[True, True, False]
)

cred_output_df.to_csv(f'conditional_chr{chr_input}_credible_sets.tsv', sep='\t', index=False)

#----------------------------------------------------------
# Iterate over SNP Data (Full)
#----------------------------------------------------------

# Create ouptut file with header from one of the input files
file_path = glob.glob('full_cred_sets/*.snp')[0]
snps = pd.read_csv(file_path, sep=' ', index_col=0)

with open(f'full_chr{chr_input}_pips.tsv', 'w') as f_out:
    header = snps.columns.tolist()
    header = ['gene'] + header + ['\n']
    f_out.write('\t'.join(header))

for file_path in sorted(glob.glob('full_cred_sets/*.snp')):

    gene_name = re.sub('\\.snp', '', os.path.basename(file_path))

    snps = pd.read_csv(file_path, sep=' ', index_col=0).sort_values(by='position')
    snps.insert(0, 'Gene', [gene_name] * snps.shape[0])

    snps.to_csv(f'full_chr{chr_input}_pips.tsv', sep='\t', header=False, index=False, mode='a')

#----------------------------------------------------------
# Iterate over SNP Data (Conditional)
#----------------------------------------------------------

# Create ouptut file with header from one of the input files
file_path = glob.glob('conditional_cred_sets/*.snp')[0]
snps = pd.read_csv(file_path, sep=' ', index_col=0)

with open(f'conditional_chr{chr_input}_pips.tsv', 'w') as f_out:
    header = snps.columns.tolist()
    header = ['gene', 'signal'] + header + ['\n']
    f_out.write('\t'.join(header))

for file_path in sorted(glob.glob('conditional_cred_sets/*.snp')):

    locus_name = re.sub('\\.snp', '', os.path.basename(file_path))
    gene_name = re.sub('-.*', '', locus_name)
    signal = re.sub('.*-', '', locus_name)

    snps = pd.read_csv(file_path, sep=' ', index_col=0).sort_values(by='position')
    snps.insert(0, 'Gene', [gene_name] * snps.shape[0])
    snps.insert(1, 'Signal', [signal] * snps.shape[0])

    snps.to_csv(f'conditional_chr{chr_input}_pips.tsv', sep='\t', header=False, index=False, mode='a')
