#----------------------------------------------------------
# Output cis-eQTL
# Created: 22 February 2022
#----------------------------------------------------------

#----------------------------------------------------------
# Environment
#----------------------------------------------------------

rm(list=ls())

library(tidyverse)
library(data.table)

options(stringsAsFactors = FALSE)

#----------------------------------------------------------
# Load Data
#----------------------------------------------------------

# Initial cis-eQTL Lead SNPs
cis.eqtl <- read.table("/nfs/users/nfs_n/nm18/gains_team282/eqtl/cisresults/eigenMT/ciseqtl_eigenMT_corrected.txt") %>%
    dplyr::filter(Sig)

# Conditional Analysis Results
cis.eqtl.conditional <- readRDS("/nfs/users/nfs_n/nm18/gains_team282/eqtl/cisresults/conditionalanalysis/conditional_eQTL_results_final.rds")

# cis-pQTL
cis.pqtl <- readRDS("/nfs/users/nfs_n/nm18/gains_team282/proteomics/pqtl/pqtl_ms2019/cis_pqtl_all.RDS") %>%
    dplyr::group_by(Gene) %>%
    dplyr::slice_min(order_by=pQTL_pval, n=1)

# trans-pQTL
trans.pqtl <- readRDS("/nfs/users/nfs_n/nm18/gains_team282/proteomics/pqtl/pqtl_ms2019/trans_pqtl_all.RDS") %>%
    dplyr::group_by(pQTL_Protein, pQTL_Locus) %>%
    dplyr::slice_min(order_by=pQTL_pval, n=1)

# module QTL
module.ss.dir = "/nfs/users/nfs_n/nm18/gains_team282/nikhil/expression/eigengene_sva/wgcna_summary_statistics/"
module.qtl <- do.call(rbind, lapply(list.files(module.ss.dir, pattern="ME_[0-9]+_[0-9]+-.*\\.tsv"), function(file.name) {

    fread(paste0(module.ss.dir, file.name)) %>%
    as.data.frame() %>%
    dplyr::select(snp=1, p=5) %>%
    dplyr::mutate(module.qtl=gsub("\\.tsv", "", file.name))
})) %>%
    dplyr::group_by(module.qtl) %>%
    dplyr::slice_min(order_by=p, n=1)

# Write SNPs
write.table(gsub("\\.", ":", unique(cis.eqtl$snps)), "lead_snps.txt", quote=F, col.names=F, row.names=F)
write.table(gsub("\\.", ":", unique(cis.eqtl.conditional$SNP)), "conditional_snps.txt", quote=F, col.names=F, row.names=F)
write.table(gsub("\\.", ":", unique(cis.pqtl$SNP)), "cis_pqtl_snps.txt", quote=F, col.names=F, row.names=F)
write.table(gsub("\\.", ":", unique(trans.pqtl$SNP)), "trans_pqtl_snps.txt", quote=F, col.names=F, row.names=F)
write.table(gsub("\\.", ":", unique(module.qtl$snp)), "module_qtl_snps.txt", quote=F, col.names=F, row.names=F)
