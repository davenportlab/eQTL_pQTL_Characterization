#!/usr/bin/env bash

export NXF_WORK=/nfs/users/nfs_n/nm18/gains_team282/epigenetics/enrichment/heritability/nextflow_work/

mkdir .nextflow.heritability/
cd .nextflow.heritability/

bsub \
    -q long \
    -o heritability_output.txt \
    -e heritability_error.txt \
    -R"select[mem>32768] rusage[mem=32768]" -M32768 \
    "nextflow ../heritability.nf"
