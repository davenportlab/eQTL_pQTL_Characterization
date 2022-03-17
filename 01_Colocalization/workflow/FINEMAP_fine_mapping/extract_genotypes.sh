#!/usr/bin/env bash

export NXF_WORK=/nfs/users/nfs_n/nm18/gains_team282/nikhil/data/nextflow_work/

mkdir .nextflow.extract_genotypes/
cd .nextflow.extract_genotypes/

bsub \
    -q long \
    -o extract_genotypes_output.txt \
    -e extract_genotypes_error.txt \
    -R"select[mem>32000] rusage[mem=32000]" -M32000 \
    "nextflow ../extract_genotypes.nf"
