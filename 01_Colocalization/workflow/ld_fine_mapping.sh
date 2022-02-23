#!/usr/bin/env bash

export NXF_WORK=~/gains_team282/nikhil/colocalization/cis_eqtl/fine_mapping/LD/nextflow_work/

mkdir .nextflow.ld_fine_mapping/
cd .nextflow.ld_fine_mapping/

bsub \
    -q long \
    -o ld_fine_mapping_output.txt \
    -e ld_fine_mapping_error.txt \
    -R"select[mem>65536] rusage[mem=65536]" -M65536 \
    "nextflow ../ld_fine_mapping.nf"

cd ../