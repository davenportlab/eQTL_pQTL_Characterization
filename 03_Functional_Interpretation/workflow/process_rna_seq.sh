#!/usr/bin/env bash

export NXF_WORK=~/gains_team282/epigenetics/calderon_et_al/processed/rna_seq/nextflow_work/

bsub -q normal -o calderon_et_al_rna_seq_output.txt -e calderon_et_al_rna_seq_error.txt -R"select[mem>8192] rusage[mem=8192]" -M8192 \
    "nextflow process_rna_seq.nf --reads_dir ~/gains_team282/epigenetics/calderon_et_al/raw/rna_seq/ --output ~/gains_team282/epigenetics/calderon_et_al/processed/rna_seq/"
