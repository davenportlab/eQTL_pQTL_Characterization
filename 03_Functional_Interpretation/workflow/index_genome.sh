#!/usr/bin/env bash

bsub -q normal -o index_genome_output.txt -e index_genome_error.txt -R"select[mem>8192] rusage[mem=8192]" -M8192 \
    "nextflow index_genome.nf"
