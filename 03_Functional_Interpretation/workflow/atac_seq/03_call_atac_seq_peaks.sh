#!/usr/bin/env bash

export NXF_WORK=/nfs/users/nfs_n/nm18/gains_team282/epigenetics/accessibility/merged/atac_seq/nextflow_work/

mkdir .nextflow.03_call_atac_seq_peaks/
cd .nextflow.03_call_atac_seq_peaks/

bsub \
    -q long \
    -o calderon_et_al_call_atac_seq_peaks_output.txt \
    -e calderon_et_al_call_atac_seq_peaks_error.txt \
    -R"select[mem>8192] rusage[mem=8192]" -M8192 \
    "nextflow ../03_call_atac_seq_peaks.nf --atac_seq_dir /nfs/users/nfs_n/nm18/gains_team282/epigenetics/accessibility/processed/atac_seq/ --output /nfs/users/nfs_n/nm18/gains_team282/epigenetics/accessibility/merged/atac_seq/"
