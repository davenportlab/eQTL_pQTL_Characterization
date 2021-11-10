#!/usr/bin/env bash

export NXF_WORK=~/gains_team282/nikhil/colocalization/cis_eqtl/nextflow_work/

bsub -q normal -o prepare_eqtl_data_for_coloc_output.txt -e prepare_eqtl_data_for_coloc_error.txt -R"select[mem>32000] rusage[mem=32000]" -M32000 \
    "nextflow prepare_eqtl_data_for_coloc.nf"
