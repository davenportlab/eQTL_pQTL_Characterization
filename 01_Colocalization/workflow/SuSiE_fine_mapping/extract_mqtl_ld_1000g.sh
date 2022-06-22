#!/usr/bin/env bash

export NXF_WORK=/nfs/users/nfs_n/nm18/gains_team282/nikhil/data/1000G/nextflow_work/

mkdir .nextflow.extract_mqtl_ld_1000g/
cd .nextflow.extract_mqtl_ld_1000g\

bsub \
    -q long \
    -o extract_mqtl_ld_1000g_output.txt \
    -e extract_mqtl_ld_1000g_error.txt \
    -R"select[mem>65536] rusage[mem=65536]" -M65536 \
    "nextflow ../extract_mqtl_ld_1000g.nf"

cd ../
