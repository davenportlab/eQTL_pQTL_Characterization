#!/usr/bin/env bash

export NXF_WORK=~/gains_team282/epigenetics/calderon_et_al/analysis/atac_seq/nextflow_work/

mkdir .nextflow.atac_seq_motifs_ia/
cd .nextflow.atac_seq_motifs_ia/

bsub \
    -q long \
    -o calderon_et_al_atac_seq_homer_output.txt \
    -e calderon_et_al_atac_seq_homer_error.txt \
    -R"select[mem>8192] rusage[mem=8192]" -M8192 \
    "nextflow ../atac_seq_motifs.nf --da_peak_set /nfs/users/nfs_n/nm18/gains_team282/epigenetics/accessibility/analysis/atac_seq/da_peak_ia_set.csv --output_dir /nfs/users/nfs_n/nm18/gains_team282/epigenetics/accessibility/analysis/atac_seq/immune_atlas/"

cd ../

mkdir .nextflow.atac_seq_motifs_rm/
cd .nextflow.atac_seq_motifs_rm/

bsub \
    -q long \
    -o calderon_et_al_atac_seq_homer_output.txt \
    -e calderon_et_al_atac_seq_homer_error.txt \
    -R"select[mem>8192] rusage[mem=8192]" -M8192 \
    "nextflow ../atac_seq_motifs.nf --da_peak_set /nfs/users/nfs_n/nm18/gains_team282/epigenetics/accessibility/analysis/atac_seq/da_peak_rm_set.csv --output_dir /nfs/users/nfs_n/nm18/gains_team282/epigenetics/accessibility/analysis/atac_seq/neutrophil_atlas/"

