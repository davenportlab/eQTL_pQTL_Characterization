#!/usr/bin/env bash

export NXF_WORK=/nfs/users/nfs_n/nm18/gains_team282/epigenetics/accessibility/analysis/atac_seq/nextflow_work/

mkdir .nextflow.04_sample_quality/
cd .nextflow.04_sample_quality/

bsub \
    -q long \
    -o sample_quality_output.txt \
    -e sample_quality_error.txt \
    -R"select[mem>8192] rusage[mem=8192]" -M8192 \
    "nextflow ../04_sample_quality.nf"
