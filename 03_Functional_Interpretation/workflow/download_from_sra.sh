#!/usr/bin/env bash

export NXF_WORK=/nfs/users/nfs_n/nm18/gains_team282/epigenetics/accessibility/nextflow_work/

mkdir .nextflow.download_from_sra/
cd .nextflow.download_from_sra/

## Download after NCBI fixes RNA-Seq Data for Ram-Mohan et al.
# bsub \
#     -q normal \
#     -o rna_seq_output.txt \
#     -e rna_seq_error.txt \
#     -R"select[mem>8192] rusage[mem=8192]" -M8192 \
#     "nextflow ../download_from_sra.nf --section accessibility --assay rna_seq --sra_table /nfs/users/nfs_n/nm18/eQTL_pQTL_Characterization/03_Functional_Interpretation/metadata/reads_rna_seq.txt"

bsub \
    -q normal \
    -o atac_seq_output.txt \
    -e atac_seq_error.txt \
    -R"select[mem>8192] rusage[mem=8192]" -M8192 \
    "nextflow ../download_from_sra.nf --section accessibility --assay atac_seq --sra_table /nfs/users/nfs_n/nm18/eQTL_pQTL_Characterization/03_Functional_Interpretation/metadata/reads_atac_seq.txt"
