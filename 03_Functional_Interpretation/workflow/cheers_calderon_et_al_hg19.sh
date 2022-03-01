#!/usr/bin/env bash

export NXF_WORK=/nfs/users/nfs_n/nm18/gains_team282/epigenetics/calderon_et_al_hg19/nextflow_work/

mkdir .nextflow.cheers_calderon_et_al_hg19/
cd .nextflow.cheers_calderon_et_al_hg19/

bsub \
    -q long \
    -o cheers_output.txt \
    -e cheers_errors.txt \
    -R"select[mem>32768] rusage[mem=32768]" -M32768 \
    "nextflow ../cheers_calderon_et_al_hg19.nf"
