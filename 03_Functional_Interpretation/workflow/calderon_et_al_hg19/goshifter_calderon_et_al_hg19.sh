#!/usr/bin/env bash

export NXF_WORK=/nfs/users/nfs_n/nm18/gains_team282/epigenetics/calderon_et_al_hg19/nextflow_work/

mkdir .nextflow.goshifter_calderon_et_al_hg19/
cd .nextflow.goshifter_calderon_et_al_hg19/

bsub \
    -q long \
    -o goshifter_output.txt \
    -e goshifter_errors.txt \
    -R"select[mem>32768] rusage[mem=32768]" -M32768 \
    "nextflow ../goshifter_calderon_et_al_hg19.nf"
