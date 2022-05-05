#!/usr/bin/env bash

export NXF_WORK=/nfs/users/nfs_n/nm18/gains_team282/nikhil/expression/eigengene_sva/nextflow_work/

mkdir .nextflow.module_qtl_summary_statistics/
cd .nextflow.module_qtl_summary_statistics/

bsub \
    -q long \
    -o module_qtl_summary_statistics_output.txt \
    -e module_qtl_summary_statistics_error.txt \
    -R"select[mem>8000] rusage[mem=8000]" -M8000 \
    "nextflow ../module_qtl_summary_statistics.nf --output_dir /nfs/users/nfs_n/nm18/gains_team282/nikhil/expression/eigengene_sva/wgcna_summary_statistics/"
