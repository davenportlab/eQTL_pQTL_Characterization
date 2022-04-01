#!/usr/bin/env bash

export NXF_WORK=/nfs/users/nfs_n/nm18/gains_team282/nikhil/colocalization/pqtl/conditional_effects/FINEMAP/nextflow_work/

mkdir .nextflow.pqtl_conditional_effects_finemap/
cd .nextflow.pqtl_conditional_effects_finemap/

bsub \
    -q long \
    -o pqtl_conditional_effects_finemap_output.txt \
    -e pqtl_conditional_effects_finemap_error.txt \
    -R"select[mem>32768] rusage[mem=32768]" -M32768 \
    "nextflow ../pqtl_conditional_effects_finemap.nf"
