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

options(stringsAsFactors = FALSE)

setwd("~/eQTL_pQTL_Characterization/")

source("01_Colocalization/scripts/utils/ggplot_theme.R")

#----------------------------------------------------------
# Load Gene Expression Data
#----------------------------------------------------------

gene.exp <- read.table("/nfs/team282/data/gains_team282/logcpm_864_20397_hla.txt")

gene.exp.sd <- data.frame(
  Gene=rownames(gene.exp),
  SD=apply(gene.exp, 1, sd)
)

saveRDS(gene.exp.sd, "~/gains_team282/nikhil/colocalization/cis_eqtl/cis.eqtl.sdY.RDS")
