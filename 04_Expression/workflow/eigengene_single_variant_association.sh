#!/usr/bin/env bash

export NXF_WORK=~/gains_team282/nikhil/expression/eigengene_sva/nextflow_work/

for ME_NUM in 87 91 100
do
    mkdir .nextflow.ME_${ME_NUM}/
    cd .nextflow.ME_${ME_NUM}/

    bsub \
        -q long \
        -o eigengene_single_variant_association_ME_${ME_NUM}_output.txt \
        -e eigengene_single_variant_association_ME_${ME_NUM}_error.txt \
        -R"select[mem>8000] rusage[mem=8000]" -M8000 \
        "nextflow ../eigengene_single_variant_association.nf --me ME_${ME_NUM}"
    
    cd ../
done
