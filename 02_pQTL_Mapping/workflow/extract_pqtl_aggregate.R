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

wg.pqtl <- do.call(rbind, mclapply(list.files("trans_pqtl/"), function(file.name) {
    readRDS(paste0("trans_pqtl/", file.name))
}, mc.cores=16))

saveRDS(wg.pqtl, "trans_pqtl_all.RDS")

cis.pqtl <- do.call(rbind, mclapply(list.files("cis_pqtl/"), function(file.name) {
    readRDS(paste0("cis_pqtl/", file.name))
}, mc.cores=16))

saveRDS(cis.pqtl, "cis_pqtl_all.RDS")
