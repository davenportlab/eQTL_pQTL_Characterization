#!/usr/bin/env bash

export NXF_WORK=/nfs/users/nfs_n/nm18/gains_team282/nikhil/colocalization/mqtl/fine_mapping/FINEMAP/nextflow_work/


mkdir .nextflow.mqtl_finemap_fine_mapping/
cd .nextflow.mqtl_finemap_fine_mapping/

bsub \
    -q long \
    -o finemap_fine_mapping_output.txt \
    -e finemap_fine_mapping_error.txt \
    -R"select[mem>65536] rusage[mem=65536]" -M65536 \
    "nextflow ../mqtl_finemap_fine_mapping.nf"
