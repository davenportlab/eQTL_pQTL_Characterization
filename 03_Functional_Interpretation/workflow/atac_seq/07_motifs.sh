#!/usr/bin/env bash

export NXF_WORK=~/gains_team282/epigenetics/accessibility/analysis/atac_seq/nextflow_work/

mkdir .nextflow.07_motifs/
cd .nextflow.07_motifs/

bsub \
    -q long \
    -o motifs_output.txt \
    -e motifs_error.txt \
    -R"select[mem>8192] rusage[mem=8192]" -M8192 \
    "nextflow ../07_motifs.nf"
