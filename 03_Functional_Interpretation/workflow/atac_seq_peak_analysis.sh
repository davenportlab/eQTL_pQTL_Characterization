#!/usr/bin/env bash

for atlas in combined immune neutrophil
do

    export NXF_WORK=/nfs/users/nfs_n/nm18/gains_team282/epigenetics/accessibility/analysis/atac_seq/$atlas/nextflow_work/

    mkdir .nextflow.atac_seq_${atlas}_peak_analysis/
    cd .nextflow.atac_seq_${atlas}_peak_analysis/

    bsub \
        -q long \
        -o calderon_et_al_atac_seq_${atlas}_peak_analysis_output.txt \
        -e calderon_et_al_atac_seq_${atlas}_peak_analysis_error.txt \
        -R"select[mem>8192] rusage[mem=8192]" -M8192 \
        "nextflow ../atac_seq_peak_analysis.nf --atlas $atlas --peak_saf /nfs/users/nfs_n/nm18/gains_team282/epigenetics/accessibility/analysis/atac_seq/$atlas/consensus_peaks.saf --output_dir /nfs/users/nfs_n/nm18/gains_team282/epigenetics/accessibility/analysis/atac_seq/$atlas/"

    cd ../

done
