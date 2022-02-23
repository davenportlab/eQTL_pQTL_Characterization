#!/usr/bin/env bash

export NXF_WORK=/nfs/users/nfs_n/nm18/gains_team282/proteomics/pqtl/pqtl_ms2019/nextflow_work/

mkdir .nextflow.extract_pqtl/
cd .nextflow.extract_pqtl/

bsub \
    -q long \
    -o extract_pqtl_output.txt \
    -e extract_pqtl_error.txt \
    -R"select[mem>65536] rusage[mem=65536]" -M65536 \
    "nextflow ../extract_pqtl.nf"