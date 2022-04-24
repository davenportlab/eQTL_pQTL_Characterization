#!/usr/bin/env bash

export NXF_WORK=/nfs/users/nfs_n/nm18/gains_team282/epigenetics/accessibility/analysis/atac_seq/nextflow_work/

mkdir .nextflow.05_consensus_peak_set/
cd .nextflow.05_consensus_peak_set/

bsub \
    -q long \
    -o consensus_peak_set_output.txt \
    -e consensus_peak_set_error.txt \
    -R"select[mem>8192] rusage[mem=8192]" -M8192 \
    "nextflow ../05_consensus_peak_set.nf"
