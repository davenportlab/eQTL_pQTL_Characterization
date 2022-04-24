#!/usr/bin/env bash

for atlas in immune neutrophil
do

    export NXF_WORK=/nfs/users/nfs_n/nm18/gains_team282/epigenetics/accessibility/analysis/atac_seq/$atlas/nextflow_work/

    mkdir .nextflow.06_count_fragments_${atlas}/
    cd .nextflow.06_count_fragments_${atlas}/

    bsub \
        -q long \
        -o count_fragments_${atlas}_output.txt \
        -e count_fragments_${atlas}_error.txt \
        -R"select[mem>8192] rusage[mem=8192]" -M8192 \
        "nextflow ../06_count_fragments.nf \\
            --atlas $atlas \\
            --peak_saf /nfs/users/nfs_n/nm18/gains_team282/epigenetics/accessibility/analysis/atac_seq/$atlas/consensus_peaks.saf \\
            --peak_set_dir /nfs/users/nfs_n/nm18/gains_team282/epigenetics/accessibility/analysis/atac_seq/$atlas/peak_sets/ \\
            --cell_type_set_dir /nfs/users/nfs_n/nm18/gains_team282/epigenetics/accessibility/analysis/atac_seq/$atlas/cell_type_peak_sets/ \\
            --output_dir /nfs/users/nfs_n/nm18/gains_team282/epigenetics/accessibility/analysis/atac_seq/$atlas/"

    cd ../

done
