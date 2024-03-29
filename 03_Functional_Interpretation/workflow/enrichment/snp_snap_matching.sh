#!/usr/bin/env bash

export NXF_WORK=/nfs/users/nfs_n/nm18/gains_team282/epigenetics/enrichment/snp_snap/nextflow_work/

mkdir .nextflow.snp_snap_matching/
cd .nextflow.snp_snap_matching/

bsub \
    -q long \
    -o snp_snap_matching_output.txt \
    -e snp_snap_matching_error.txt \
    -R"select[mem>32768] rusage[mem=32768]" -M32768 \
    "nextflow ../snp_snap_matching.nf"
