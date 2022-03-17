#!/usr/bin/env bash

export NXF_WORK=~/gains_team282/epigenetics/accessibility/analysis/atac_seq/nextflow_work/

mkdir .nextflow.06_atac_seq_motifs/
cd .nextflow.06_atac_seq_motifs/

bsub \
    -q long \
    -o calderon_et_al_atac_seq_homer_output.txt \
    -e calderon_et_al_atac_seq_homer_error.txt \
    -R"select[mem>8192] rusage[mem=8192]" -M8192 \
    "nextflow ../06_atac_seq_motifs.nf"
