#!/usr/bin/env bash

export NXF_WORK=~/gains_team282/nikhil/colocalization/cis_eqtl/nextflow_work/

bsub -q long -o prepare_eqtl_data_for_coloc_output_X.txt -e prepare_eqtl_data_for_coloc_error_X.txt -R"select[mem>65536] rusage[mem=65536]" -M65536 \
    "nextflow prepare_eqtl_data_for_coloc.nf --chr X"
