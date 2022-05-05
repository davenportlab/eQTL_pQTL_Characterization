#!/usr/bin/env bash

export NXF_WORK=~/gains_team282/nikhil/data/nextflow_work/

mkdir .nextflow.extract_ss_genotypes/
cd .nextflow.extract_ss_genotypes/

bsub \
    -q long \
    -o extract_ss_genotypes_output.txt \
    -e extract_ss_genotypes_error.txt \
    -R"select[mem>32000] rusage[mem=32000]" -M32000 \
    "nextflow ../extract_ss_genotypes.nf"
