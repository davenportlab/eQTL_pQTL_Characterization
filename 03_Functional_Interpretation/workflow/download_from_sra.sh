#!/usr/bin/env bash

export NXF_WORK=~/gains_team282/epigenetics/nextflow_work/

bsub -q normal -o calderon_et_al_rna_seq_output.txt -e calderon_et_al_rna_seq_error.txt -R"select[mem>8192] rusage[mem=8192]" -M8192 \
    "nextflow download_from_sra.nf --study calderon_et_al --assay rna_seq --library_layout paired --sra_table ../metadata/reads_calderon_et_al_rna_seq.txt"

bsub -q normal -o calderon_et_al_chip_seq_output.txt -e calderon_et_al_chip_seq_error.txt -R"select[mem>8192] rusage[mem=8192]" -M8192 \
    "nextflow download_from_sra.nf --study calderon_et_al --assay chip_seq --library_layout paired --sra_table ../metadata/reads_calderon_et_al_chip_seq.txt"

bsub -q normal -o calderon_et_al_atac_seq_output.txt -e calderon_et_al_atac_seq_error.txt -R"select[mem>8192] rusage[mem=8192]" -M8192 \
    "nextflow download_from_sra.nf --study calderon_et_al --assay atac_seq --library_layout paired --sra_table ../metadata/reads_calderon_et_al_atac_seq.txt"

bsub -q normal -o brands_et_al_rna_seq_output.txt -e brands_et_al_rna_seq_error.txt -R"select[mem>8192] rusage[mem=8192]" -M8192 \
    "nextflow download_from_sra.nf --study brands_et_al --assay rna_seq --library_layout single --sra_table ../metadata/reads_brands_et_al_rna_seq.txt"

bsub -q normal -o brands_et_al_bisulfite_seq_output.txt -e brands_et_al_bisulfite_seq_error.txt -R"select[mem>8192] rusage[mem=8192]" -M8192 \
    "nextflow download_from_sra.nf --study brands_et_al --assay bisulfite_seq --library_layout single --sra_table ../metadata/reads_brands_et_al_bisulfite_seq.txt"
