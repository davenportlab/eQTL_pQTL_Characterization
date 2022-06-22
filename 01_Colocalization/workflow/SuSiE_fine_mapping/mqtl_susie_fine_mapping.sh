#!/usr/bin/env bash

export NXF_WORK=/nfs/users/nfs_n/nm18/gains_team282/nikhil/colocalization/mqtl/fine_mapping/SuSiE/nextflow_work/

mkdir .nextflow.mqtl_susie_fine_mapping/
cd .nextflow.mqtl_susie_fine_mapping\

bsub \
    -q long \
    -o susie_fine_mapping_output.txt \
    -e susie_fine_mapping_error.txt \
    -R"select[mem>65536] rusage[mem=65536]" -M65536 \
    "nextflow ../mqtl_susie_fine_mapping.nf"

cd ../
