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
library(parallel)
library(foreach)

options(stringsAsFactors = FALSE)

#----------------------------------------------------------
# Load Data
#----------------------------------------------------------

# Initial cis-eQTL Lead SNPs
cis.eqtl <- read.table("/nfs/users/nfs_n/nm18/gains_team282/eqtl/cisresults/eigenMT/ciseqtl_eigenMT_corrected.txt") %>%
    dplyr::filter(Sig)

# Conditional Analysis Results
cis.eqtl.conditional <- readRDS("/nfs/users/nfs_n/nm18/gains_team282/eqtl/cisresults/conditionalanalysis/conditional_eQTL_results_final.rds")

# Write SNPs
write.table(unique(cis.eqtl$snps), "lead_snps.txt", quote=F, col.names=F, row.names=F)
write.table(unique(cis.eqtl.conditional$SNP), "conditional_snps.txt", quote=F, col.names=F, row.names=F)
