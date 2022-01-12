#!/usr/bin/env bash

export NXF_WORK=~/gains_team282/nikhil/expression/eigengene_sva/nextflow_work/

for ME_NUM in {52..102}
do
    mkdir .nextflow.ME_${ME_NUM}/
    cd .nextflow.ME_${ME_NUM}/

    bsub \
        -q long \
        -o eigengene_single_variant_association_ME_${ME_NUM}_output.txt \
        -e eigengene_single_variant_association_ME_${ME_NUM}_error.txt \
        -R"select[mem>8000] rusage[mem=8000]" -M8000 \
        "nextflow ../eigengene_single_variant_association.nf --me ME_${ME_NUM} --design_matrix /nfs/users/nfs_n/nm18/gains_team282/nikhil/expression/eigengene_sva/mapping_data_cell_proportions.csv --output_dir /nfs/users/nfs_n/nm18/gains_team282/nikhil/expression/eigengene_sva/cell_proportions/"
    
    cd ../
done
