#!/usr/bin/env bash

export NXF_WORK=/nfs/users/nfs_n/nm18/gains_team282/nikhil/colocalization/cis_eqtl/fine_mapping/FINEMAP/nextflow_work/

for CHR in 1 2 3 4 6
do

    mkdir .nextflow.cis_eqtl_finemap_fine_mapping_CHR_${CHR}/
    cd .nextflow.cis_eqtl_finemap_fine_mapping_CHR_${CHR}/

    bsub \
        -q long \
        -o finemap_fine_mapping_output_${CHR}.txt \
        -e finemap_fine_mapping_error_${CHR}.txt \
        -R"select[mem>65536] rusage[mem=65536]" -M65536 \
        "nextflow ../cis_eqtl_finemap_fine_mapping.nf --chr $CHR"

    cd ../

done
