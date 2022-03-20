#!/usr/bin/env bash

export NXF_WORK=/nfs/users/nfs_n/nm18/gains_team282/epigenetics/accessibility/analysis/atac_seq/nextflow_work/

mkdir .nextflow.04_atac_seq_sample_peak_analysis/
cd .nextflow.04_atac_seq_sample_peak_analysis/

bsub \
    -q long \
    -o atac_seq_sample_peak_analysis_output.txt \
    -e atac_seq_sample_peak_analysis_error.txt \
    -R"select[mem>8192] rusage[mem=8192]" -M8192 \
    "nextflow ../04_atac_seq_sample_peak_analysis.nf"
