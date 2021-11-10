#----------------------------------------------------------
# COLOC with SuSiE
# Created: 25 October 2021
#----------------------------------------------------------

#----------------------------------------------------------
# Environment
#----------------------------------------------------------

rm(list=ls())

library(tidyverse)
library(coloc)
library(susieR)

options(stringsAsFactors = FALSE)

setwd("~/eQTL_pQTL_Characterization/")

source("01_Colocalization/scripts/utils/ggplot_theme.R")

#----------------------------------------------------------
# Load Data
#----------------------------------------------------------

cis.eqtl.loci <- readRDS("~/gains_team282/nikhil/colocalization/cis_eqtl/cis.eqtl.loci.RDS")

cis.pqtl.loci <- readRDS("~/gains_team282/nikhil/colocalization/pilot.cis.pqtl.loci.RDS")
cis.pqtl.sdY <- readRDS("~/gains_team282/nikhil/colocalization/pilot.cis.pqtl.sdY.RDS")

#----------------------------------------------------------
# Run Colocalization Analysis
#----------------------------------------------------------

pqtl.ensembl.genes <- lapply(cis.pqtl.loci, function(x) strsplit(x$Ensembl.ID[1], "; "))
eqtl.to.test <- lapply(pqtl.ensembl.genes, function(x) names(cis.eqtl.loci)[names(cis.eqtl.loci) %in% x])

results <- list()
results[["pQTL Locus"]] = list()
results[["eQTL Locus"]] = list()
results[["PP H3"]] = list()
results[["PP H4"]] = list()
counter <- 0

for (pqtl in names(eqtl.to.test)) {
  
  for (eqtl in eqtl.to.test[[pqtl]]) {
    
    counter <- counter + 1
    
    cis.eqtl.locus = cis.eqtl.loci[[eqtl]]
    cis.pqtl.locus = cis.pqtl.loci[[pqtl]]
    
    rownames(cis.pqtl.locus) = cis.pqtl.locus$snp

    snps.common = intersect(cis.eqtl.locus$snp, cis.pqtl.locus$snp)

    filter = cis.eqtl.locus$snp %in% snps.common & !apply(cis.eqtl.locus$LD, 1, function(x) length(which(x < 0)) == dim(cis.eqtl.locus$LD)[1] - 1)
    cis.eqtl.locus$snp = cis.eqtl.locus$snp[filter]
    cis.eqtl.locus$position = cis.eqtl.locus$position[filter]
    cis.eqtl.locus$beta = cis.eqtl.locus$beta[filter]
    cis.eqtl.locus$varbeta = cis.eqtl.locus$varbeta[filter]
    cis.eqtl.locus$LD = cis.eqtl.locus$LD[filter, filter]
    
    cis.pqtl.locus = cis.pqtl.locus[cis.eqtl.locus$snp,] %>% as.list()
    cis.pqtl.locus$position = cis.eqtl.locus$position
    cis.pqtl.locus$type = "quant"
    cis.pqtl.locus$sdY = cis.pqtl.sdY[pqtl, "SD"]
    cis.pqtl.locus$LD = cis.eqtl.locus$LD

    if (length(cis.pqtl.locus$snp) == 0) {
      next
    }
    
    result <- coloc.susie(cis.eqtl.locus, cis.pqtl.locus)
    
    results[["pQTL Locus"]][[counter]] <- pqtl
    results[["eQTL Locus"]][[counter]] <- eqtl
    results[["Summary"]][[counter]] <- result$summary
  }
}

saveRDS("~/gains_team282/nikhil/colocalization/cis.eqtl.pilot.cis.pqtl.colocalization.RDS")
