#!/usr/bin/env bash

export NXF_WORK=/nfs/users/nfs_n/nm18/gains_team282/nikhil/colocalization/cis_eqtl/conditional_effects/LD/nextflow_work/

for CHR in {1..6}
do

    mkdir .nextflow.eqtl_conditional_effects_ld_chr_${CHR}/
    cd .nextflow.eqtl_conditional_effects_ld_chr_${CHR}/

    bsub \
        -q long \
        -o eqtl_conditional_effects_ld_output.txt \
        -e eqtl_conditional_effects_ld_error.txt \
        -R"select[mem>32768] rusage[mem=32768]" -M32768 \
        "nextflow ../eqtl_conditional_effects_ld.nf --chr $CHR"

    cd ../

done
