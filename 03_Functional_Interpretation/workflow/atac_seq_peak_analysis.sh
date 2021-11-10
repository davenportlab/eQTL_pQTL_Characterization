#!/usr/bin/env bash

export NXF_WORK=~/gains_team282/epigenetics/calderon_et_al/analysis/atac_seq/nextflow_work/

bsub -q normal -o calderon_et_al_atac_seq_peak_analysis_output.txt -e calderon_et_al_atac_seq_peak_analysis_error.txt -R"select[mem>8192] rusage[mem=8192]" -M8192 \
    "nextflow atac_seq_peak_analysis.nf --atac_seq_dir /nfs/users/nfs_n/nm18/gains_team282/epigenetics/calderon_et_al/processed/atac_seq/ --output_dir /nfs/users/nfs_n/nm18/gains_team282/epigenetics/calderon_et_al/analysis/atac_seq/"
