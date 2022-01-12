#!/usr/bin/env bash

export NXF_WORK=~/gains_team282/nikhil/colocalization/cis_eqtl/nextflow_work/

for CHR in 2 3 4
do

    bsub -q long -o susie_fine_mapping_output_${CHR}.txt -e susie_fine_mapping_error_${CHR}.txt -R"select[mem>65536] rusage[mem=65536]" -M65536 \
        "nextflow susie_fine_mapping.nf --chr $CHR"

done
