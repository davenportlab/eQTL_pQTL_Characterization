#!/usr/bin/env bash

export NXF_WORK=/nfs/users/nfs_n/nm18/gains_team282/epigenetics/enrichment/cheers/nextflow_work/

mkdir .nextflow.cheers/
cd .nextflow.cheers/

bsub \
    -q long \
    -o cheers_output.txt \
    -e cheers_errors.txt \
    -R"select[mem>32768] rusage[mem=32768]" -M32768 \
    "nextflow ../cheers.nf"
