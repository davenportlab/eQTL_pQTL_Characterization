#!/usr/bin/env bash

export NXF_WORK=~/gains_team282/nikhil/expression/eigengene_sva/nextflow_work/

bsub -q long -o extract_genotypes_output.txt -e extract_genotypes_error.txt -R"select[mem>32000] rusage[mem=32000]" -M32000 \
    "nextflow extract_genotypes.nf"
