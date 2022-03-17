#!/usr/bin/env bash

export NXF_WORK=/nfs/users/nfs_n/nm18/gains_team282/epigenetics/accessibility/merged/rna_seq/nextflow_work/

mkdir .nextflow.03_merge_replicates/
cd .nextflow.03_merge_replicates/

bsub \
    -q long \
    -o calderon_et_al_merge_replicates_output.txt \
    -e calderon_et_al_merge_replicates_error.txt \
    -R"select[mem>8192] rusage[mem=8192]" -M8192 \
    "nextflow ../03_merge_replicates.nf --rna_seq_dir /nfs/users/nfs_n/nm18/gains_team282/epigenetics/accessibility/processed/rna_seq/ --output /nfs/users/nfs_n/nm18/gains_team282/epigenetics/accessibility/merged/rna_seq/"
