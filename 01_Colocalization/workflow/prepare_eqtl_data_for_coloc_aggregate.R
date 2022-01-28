#----------------------------------------------------------
# Prepare GAinS eQTL for Colocalization
# Created: 25 October 2021
#----------------------------------------------------------

#----------------------------------------------------------
# Environment
#----------------------------------------------------------

rm(list=ls())

library(tidyverse)
library(data.table)
library(parallel)

options(stringsAsFactors = FALSE)

args = commandArgs(trailingOnly = TRUE)
chr = args[1]

#----------------------------------------------------------
# Load LD Data
#----------------------------------------------------------

cis.eqtl.sdY <- readRDS("/nfs/users/nfs_n/nm18/gains_team282/nikhil/colocalization/cis_eqtl/cis.eqtl.sdY.RDS")
locus.data.files <- list.files("locus_data/", pattern="*.csv")

cis.eqtl.loci <- mclapply(locus.data.files, function(locus.data.file) {

  locus.name = strsplit(locus.data.file, "\\.")[[1]][1]

  # Read locus data (SNP list) and locus LD matrix
  locus.data = fread(paste0("locus_data/", locus.name, ".csv"), sep=",")
  locus.ld = fread(paste0("locus_ld/", locus.name, "_ld.tsv"), sep="\t", header=TRUE)[,-1]

  # SNPs may be pruned when LD matrix is calculated
  #   Ensure that the SNP list reflects this pruning
  #   Also arrange SNPs by ascending order of position
  cis.eqtl.locus = as.data.frame(locus.data) %>%
    dplyr::filter(as.character(position) %in% colnames(locus.ld)) %>%
    dplyr::arrange(position)

  # Reorder LD matrix by ascending order of position to match SNP list
  locus.ld.order = order(as.numeric(colnames(locus.ld)))
  locus.ld.mtx = as.matrix(locus.ld)
  rownames(locus.ld.mtx) <- colnames(locus.ld.mtx)
  locus.ld.mtx = locus.ld.mtx[locus.ld.order, locus.ld.order]

  # Assign SNP names to the LD matrix
  rownames(locus.ld.mtx) <- cis.eqtl.locus$snp
  colnames(locus.ld.mtx) <- cis.eqtl.locus$snp

  # Prepare data object for COLOC
  cis.eqtl.locus = as.list(cis.eqtl.locus)[c("snp", "position", "beta", "varbeta")]
    
  cis.eqtl.locus$type = "quant"
  cis.eqtl.locus$sdY = cis.eqtl.sdY[locus.name, "SD"]

  cis.eqtl.locus$LD = locus.ld.mtx

  return(cis.eqtl.locus)
})

names(cis.eqtl.loci) <- sapply(strsplit(locus.data.files, "\\."), function(x) x[1])

saveRDS(cis.eqtl.loci, paste0("/nfs/users/nfs_n/nm18/gains_team282/nikhil/colocalization/cis_eqtl/cis.eqtl.loci.chr", chr, ".RDS"))
