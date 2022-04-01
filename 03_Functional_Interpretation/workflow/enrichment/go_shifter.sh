#!/usr/bin/env bash

export NXF_WORK=/nfs/users/nfs_n/nm18/gains_team282/epigenetics/enrichment/go_shifter/nextflow_work/

mkdir .nextflow.go_shifter/
cd .nextflow.go_shifter/

bsub \
    -q long \
    -o go_shifter_output.txt \
    -e go_shifter_errors.txt \
    -R"select[mem>32768] rusage[mem=32768]" -M32768 \
    "nextflow ../go_shifter.nf"
