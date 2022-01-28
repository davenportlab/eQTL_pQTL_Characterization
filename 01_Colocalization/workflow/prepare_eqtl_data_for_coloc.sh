#!/usr/bin/env bash

export NXF_WORK=/nfs/users/nfs_n/nm18/gains_team282/nikhil/colocalization/cis_eqtl/nextflow_work/

for CHR in 1
do

    mkdir .nextflow.CHR_${CHR}/
    cd .nextflow.CHR_${CHR}/

    bsub \
        -q long \
        -o prepare_eqtl_data_for_coloc_output_${CHR}.txt \
        -e prepare_eqtl_data_for_coloc_error_${CHR}.txt \
        -R"select[mem>65536] rusage[mem=65536]" -M65536 \
        "nextflow ../prepare_eqtl_data_for_coloc.nf --chr $CHR"

    cd ../

done
