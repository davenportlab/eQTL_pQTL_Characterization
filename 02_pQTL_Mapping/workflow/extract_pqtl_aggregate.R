#----------------------------------------------------------
# Aggregate pQTL from Summary Results
# Created: 19 February 2022
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

wg.pqtl <- do.call(rbind, mclapply(list.files("wg_pqtl/"), function(file.name) {
    readRDS(paste0("wg_pqtl/", file.name))
}, mc.cores=16))

saveRDS(wg.pqtl, "whole_genome_pqtl_all.RDS")

cis.pqtl <- do.call(rbind, mclapply(list.files("cis_pqtl/"), function(file.name) {
    readRDS(paste0("cis_pqtl/", file.name))
}, mc.cores=16))

saveRDS(cis.pqtl, "cis_pqtl_all.RDS")
