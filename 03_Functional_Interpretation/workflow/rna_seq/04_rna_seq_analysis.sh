#!/usr/bin/env bash

export NXF_WORK=/nfs/users/nfs_n/nm18/gains_team282/epigenetics/accessibility/analysis/rna_seq/nextflow_work/

mkdir .nextflow.04_rna_seq_analysis/
cd .nextflow.04_rna_seq_analysis/

bsub \
    -q long \
    -o calderon_et_al_rna_seq_analysis_output.txt \
    -e calderon_et_al_rna_seq_analysis_error.txt \
    -R"select[mem>8192] rusage[mem=8192]" -M8192 \
    "nextflow ../04_rna_seq_analysis.nf --rna_seq_dir /nfs/users/nfs_n/nm18/gains_team282/epigenetics/accessibility/merged/rna_seq/ --output /nfs/users/nfs_n/nm18/gains_team282/epigenetics/accessibility/analysis/rna_seq/"
