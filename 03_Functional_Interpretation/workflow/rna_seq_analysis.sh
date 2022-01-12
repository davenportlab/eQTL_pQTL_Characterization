#!/usr/bin/env bash

export NXF_WORK=~/gains_team282/epigenetics/calderon_et_al/analysis/rna_seq/nextflow_work/

bsub -q normal -o calderon_et_al_rna_seq_analysis_output.txt -e calderon_et_al_rna_seq_analysis_error.txt -R"select[mem>8192] rusage[mem=8192]" -M8192 \
    "nextflow rna_seq_analysis.nf --rna_seq_dir /nfs/users/nfs_n/nm18/gains_team282/epigenetics/calderon_et_al/processed/rna_seq/ --output_dir /nfs/users/nfs_n/nm18/gains_team282/epigenetics/calderon_et_al/analysis/rna_seq/"
