#----------------------------------------------------------
# Fine Mapping with SuSiE
# Created: 14 December 2021
#----------------------------------------------------------

#----------------------------------------------------------
# Environment
#----------------------------------------------------------

rm(list=ls())

library(tidyverse)
library(data.table)
library(coloc)
library(susieR)
library(parallel)
library(foreach)

options(stringsAsFactors = FALSE)

args = commandArgs(trailingOnly = TRUE)
chr.input = args[1]

#----------------------------------------------------------
# Load Data
#----------------------------------------------------------

cis.eqtl <- readRDS("/lustre/scratch119/realdata/mdt3/projects/gains_team282/eqtl/cisresults/cisqtl_all_significant.rds")

chr <- readRDS(paste0("/nfs/users/nfs_n/nm18/gains_team282/nikhil/colocalization/cis_eqtl/cis.eqtl.loci.chr", chr.input, ".RDS"))

#----------------------------------------------------------
# Split Loci in Chromosome
#----------------------------------------------------------

# Only perform fine mapping for significant cis-eQTL
loci <- intersect(names(chr), cis.eqtl$gene)

foreach(locus=loci) %dopar% {

    eqtl.locus <- chr[[locus]]
    saveRDS(eqtl.locus, paste0(locus, ".RDS"))
}
