#!/usr/bin/env bash

export NXF_WORK=/nfs/users/nfs_n/nm18/gains_team282/epigenetics/calderon_et_al/processed/atac_seq/nextflow_work/

mkdir .nextflow.process_paired_atac_seq/
cd .nextflow.process_paired_atac_seq/

bsub \
    -q long \
    -o calderon_et_al_atac_seq_processing_output.txt \
    -e calderon_et_al_atac_seq_processing_error.txt \
    -R"select[mem>8192] rusage[mem=8192]" -M8192 \
    "nextflow ../process_paired_atac_seq.nf --reads_dir /nfs/users/nfs_n/nm18/gains_team282/epigenetics/calderon_et_al/raw/atac_seq/ --output /nfs/users/nfs_n/nm18/gains_team282/epigenetics/calderon_et_al/processed/atac_seq/"
