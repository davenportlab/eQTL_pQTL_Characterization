#----------------------------------------------------------
# COLOC with One Causal Variant (OCV) Assumption
# Created: 25 October 2021
#----------------------------------------------------------

#----------------------------------------------------------
# Environment
#----------------------------------------------------------

rm(list=ls())

library(tidyverse)

if (!require("coloc")) {
  library(remotes)
  remotes::install_github("chr1swallace/coloc", upgrade=FALSE)
}

library(coloc)

options(stringsAsFactors = FALSE)

setwd("~/eQTL_pQTL_Characterization/")

source("01_Colocalization/scripts/utils/ggplot_theme.R")

#----------------------------------------------------------
# Load Data
#----------------------------------------------------------

cis.eqtl.loci <- readRDS("~/gains_team282/nikhil/colocalization/cis.eqtl.loci.RDS")
cis.eqtl.sdY <- readRDS("~/gains_team282/nikhil/colocalization/cis.eqtl.sdY.RDS")

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
    
    cis.eqtl.locus = cis.eqtl.loci[[eqtl]] %>%
      as.data.frame()
    cis.pqtl.locus = cis.pqtl.loci[[pqtl]]
    
    rownames(cis.eqtl.locus) = cis.eqtl.locus$snp
    rownames(cis.pqtl.locus) = cis.pqtl.locus$snp
    
    snps.common = intersect(cis.eqtl.locus$snp, cis.pqtl.locus$snp)
    
    cis.eqtl.locus = cis.eqtl.locus[snps.common,] %>% as.list()
    cis.pqtl.locus = cis.pqtl.locus[snps.common,] %>% as.list()
    
    cis.eqtl.locus$type = "quant"
    cis.eqtl.locus$sdY = cis.eqtl.sdY[eqtl, "SD"]
    
    cis.pqtl.locus$position = cis.eqtl.locus$position
    cis.pqtl.locus$type = "quant"
    cis.pqtl.locus$sdY = cis.pqtl.sdY[pqtl, "SD"]
    
    result <- coloc.abf(cis.eqtl.locus, cis.pqtl.locus)
    
    results[["pQTL Locus"]][[counter]] <- pqtl
    results[["eQTL Locus"]][[counter]] <- eqtl
    results[["PP H3"]][[counter]] <- result$summary[5]
    results[["PP H4"]][[counter]] <- result$summary[6]
  }
}

lapply(results, unlist) %>% as.data.frame %>%
  dplyr::mutate(PP.H3.and.H4 = PP.H3 + PP.H4) %>%
  dplyr::mutate(PP.Ratio = PP.H4 / PP.H3) %>%
  dplyr::filter(PP.H3.and.H4 >= 0.15) %>%
  dplyr::arrange(desc(PP.Ratio)) %>%
  head(n=20) %>% View
