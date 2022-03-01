#!/usr/bin/env bash

export NXF_WORK=/nfs/users/nfs_n/nm18/gains_team282/epigenetics/calderon_et_al_hg19/nextflow_work/

mkdir .nextflow.cheer_calderon_et_al_hg19/
cd .nextflow.cheer_calderon_et_al_hg19/

bsub \
    -q long \
    -o cheers.txt \
    -e cheers.txt \
    -R"select[mem>32768] rusage[mem=32768]" -M32768 \
    "nextflow ../cheers_calderon_et_al_hg19.nf"
