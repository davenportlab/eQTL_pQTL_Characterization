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
# Load eQTL Data
#----------------------------------------------------------

cis.eqtl <- readRDS("~/gains_team282/eqtl/cisresults/ciseqtl_all.rds")
geno.bim <- fread("~/gains_team282/Genotyping/All_genotyping_merged_filtered_b38_refiltered_rsID.bim", sep="\t")

colnames(geno.bim) <- c("chr", "snp", "cM", "Position", "minor_allele", "major_allele")

setnames(cis.eqtl, "snps", "snp")
setnames(cis.eqtl, "SNPpos", "position")

cis.eqtl <- merge(cis.eqtl, geno.bim, by="snp")

cis.eqtl[, chr := replace(chr.y, chr.y == 23, "X")]
cis.eqtl[, varbeta := se^2]

cis.eqtl.loci <- split(cis.eqtl[,c("snp", "position", "beta", "varbeta", "chr", "minor_allele", "major_allele")], cis.eqtl$gene)

saveRDS(cis.eqtl.loci, "~/gains_team282/nikhil/colocalization/cis.eqtl.loci.RDS")

#----------------------------------------------------------
# Load Gene Expression Data
#----------------------------------------------------------

gene.exp <- read.table("04_Expression/data/gene_expression/Logcpm_864_20417_HLA.txt")

gene.exp.sd <- data.frame(
  Gene=rownames(gene.exp),
  SD=apply(gene.exp, 1, sd)
)

gene.exp.sd <- gene.exp.sd[names(cis.eqtl.loci),]

saveRDS(gene.exp.sd, "~/gains_team282/nikhil/colocalization/cis.eqtl.sdY.RDS")
