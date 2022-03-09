#!/usr/bin/env bash

export NXF_WORK=/nfs/users/nfs_n/nm18/gains_team282/nikhil/expression/eigengene_sva/nextflow_work/

mkdir .nextflow.eigengene_single_variant_association/
cd .nextflow.eigengene_single_variant_association/

bsub \
    -q long \
    -o eigengene_single_variant_association_output.txt \
    -e eigengene_single_variant_association_error.txt \
    -R"select[mem>8000] rusage[mem=8000]" -M8000 \
    "nextflow ../eigengene_single_variant_association.nf"
