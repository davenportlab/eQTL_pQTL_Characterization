#!/usr/bin/env bash

for atlas in immune neutrophil
do

    export NXF_WORK=/nfs/users/nfs_n/nm18/gains_team282/epigenetics/regulation/$atlas/nextflow_work/

    mkdir .nextflow.peak_shape_features_${atlas}/
    cd .nextflow.peak_shape_features_${atlas}/

    bsub \
        -q long \
        -o peak_shape_features_${atlas}_output.txt \
        -e peak_shape_features_${atlas}_error.txt \
        -R"select[mem>8192] rusage[mem=8192]" -M8192 \
        "nextflow ../peak_shape_features.nf \\
            --atlas $atlas \\
            --peak_dir /nfs/users/nfs_n/nm18/gains_team282/epigenetics/accessibility/analysis/atac_seq/$atlas/cell_type_peak_sets/ \\
            --output_dir /nfs/users/nfs_n/nm18/gains_team282/epigenetics/regulation/$atlas/"

    cd ../

done
