#!/usr/bin/env bash

for atlas in immune neutrophil
do

    export NXF_WORK=/nfs/users/nfs_n/nm18/gains_team282/epigenetics/accessibility/analysis/atac_seq/$atlas/nextflow_work/

    mkdir .nextflow.07_motifs_${atlas}/
    cd .nextflow.07_motifs_${atlas}/

    bsub \
        -q long \
        -o motifs_${atlas}_output.txt \
        -e motifs_${atlas}_error.txt \
        -R"select[mem>8192] rusage[mem=8192]" -M8192 \
        "nextflow ../07_motifs.nf --atlas $atlas"

    cd ../

done
