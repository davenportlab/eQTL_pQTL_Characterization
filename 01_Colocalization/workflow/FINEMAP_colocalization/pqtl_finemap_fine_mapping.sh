#!/usr/bin/env bash

export NXF_WORK=~/gains_team282/nikhil/colocalization/pqtl/fine_mapping/FINEMAP/nextflow_work/


mkdir .nextflow.pqtl_finemap_fine_mapping/
cd .nextflow.pqtl_finemap_fine_mapping/

bsub \
    -q long \
    -o finemap_fine_mapping_output.txt \
    -e finemap_fine_mapping_error.txt \
    -R"select[mem>65536] rusage[mem=65536]" -M65536 \
    "nextflow ../pqtl_finemap_fine_mapping.nf"
