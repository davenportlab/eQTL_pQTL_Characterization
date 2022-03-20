#!/usr/bin/env bash

export NXF_WORK=/nfs/users/nfs_n/nm18/gains_team282/epigenetics/accessibility/processed/atac_seq/nextflow_work/

mkdir .nextflow.01_pre_process_paired_atac_seq/
cd .nextflow.01_pre_process_paired_atac_seq/

bsub \
    -q long \
    -o atac_seq_processing_output.txt \
    -e atac_seq_processing_error.txt \
    -R"select[mem>8192] rusage[mem=8192]" -M8192 \
    "nextflow ../01_pre_process_paired_atac_seq.nf --reads_dir /nfs/users/nfs_n/nm18/gains_team282/epigenetics/accessibility/raw/atac_seq/ --output /nfs/users/nfs_n/nm18/gains_team282/epigenetics/accessibility/processed/atac_seq/"
