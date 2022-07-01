#!/usr/bin/env bash

export NXF_WORK=/nfs/users/nfs_n/nm18/gains_team282/epigenetics/enrichment/partitioned_heritability/nextflow_work/

mkdir .nextflow.partitioned_heritability/
cd .nextflow.partitioned_heritability/

bsub \
    -q long \
    -o partitioned_heritability_output.txt \
    -e partitioned_heritability_error.txt \
    -R"select[mem>32768] rusage[mem=32768]" -M32768 \
    "nextflow ../partitioned_heritability.nf"
