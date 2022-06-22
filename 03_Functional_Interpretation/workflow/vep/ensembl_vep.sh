#!/usr/bin/env bash

cd /nfs/users/nfs_n/nm18/gains_team282/epigenetics/vep/ensembl_vep/

# Remove logs from previous run
rm -f ensembl_vep_output.txt ensembl_vep_error.txt

# Remove results from previous run
rm -f qtl_both_ref.vep.vcf.gz qtl_both_ref.vep.vcf.gz_summary.html

# Regenerate QTL VCF file
rm -f qtl_both_ref.vcf.gz
sed 's/^chr//g' /nfs/users/nfs_n/nm18/gains_team282/epigenetics/vep/qtl_both_ref.vcf > qtl_both_ref.vcf
gzip qtl_both_ref.vcf

# Execute HGI's VEP
bsub \
    -q normal \
    -o ensembl_vep_output.txt \
    -e ensembl_vep_error.txt \
    -R"select[mem>32768] rusage[mem=32768]" -M32768 \
    "/lustre/scratch118/humgen/resources/ensembl/vep/run_vep_104.0.sh qtl_both_ref.vcf.gz"
