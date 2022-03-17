#!/usr/bin/env bash

export NXF_WORK=/nfs/users/nfs_n/nm18/gains_team282/epigenetics/accessibility/processed/rna_seq/nextflow_work/

mkdir .nextflow.02_process_paired_rna_seq/
cd .nextflow.02_process_paired_rna_seq/

bsub \
    -q long \
    -o calderon_et_al_rna_seq_processing_output.txt \
    -e calderon_et_al_rna_seq_processing_error.txt \
    -R"select[mem>8192] rusage[mem=8192]" -M8192 \
    "nextflow ../02_process_paired_rna_seq.nf --output /nfs/users/nfs_n/nm18/gains_team282/epigenetics/accessibility/processed/rna_seq/"
