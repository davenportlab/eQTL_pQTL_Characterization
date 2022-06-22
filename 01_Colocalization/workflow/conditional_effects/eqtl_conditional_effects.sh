#!/usr/bin/env bash

export NXF_WORK=/nfs/users/nfs_n/nm18/gains_team282/nikhil/colocalization/cis_eqtl/conditional_effects/LMM/nextflow_work/

for CHR in {1..5}
do

    mkdir .nextflow.eqtl_conditional_effects_chr_${CHR}/
    cd .nextflow.eqtl_conditional_effects_chr_${CHR}/

    bsub \
        -q long \
        -o eqtl_conditional_effects_output.txt \
        -e eqtl_conditional_effects_error.txt \
        -R"select[mem>32768] rusage[mem=32768]" -M32768 \
        "nextflow ../eqtl_conditional_effects.nf --chr $CHR --genotypes /nfs/users/nfs_n/nm18/gains_team282/nikhil/data/genotypes/eqtl_genotypes_${CHR}.raw"

    cd ../

done
