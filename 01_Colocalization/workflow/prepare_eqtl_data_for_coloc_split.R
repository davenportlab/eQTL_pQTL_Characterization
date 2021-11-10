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

#----------------------------------------------------------
# Load Gene Expression Data
#----------------------------------------------------------

gene.exp <- read.table("/nfs/team282/data/gains_team282/logcpm_864_20397_hla.txt")

gene.exp.sd <- data.frame(
  Gene=rownames(gene.exp),
  SD=apply(gene.exp, 1, sd)
)

saveRDS(gene.exp.sd, "~/gains_team282/nikhil/colocalization/cis_eqtl/cis.eqtl.sdY.RDS")

#----------------------------------------------------------
# Load eQTL Data
#----------------------------------------------------------

cis.eqtl <- readRDS("~/gains_team282/eqtl/cisresults/ciseqtl_all.rds")
geno.bim <- fread("~/gains_team282/Genotyping/All_genotyping_merged_filtered_b38_refiltered_rsID.bim", sep="\t")

cat("Data Loaded into Memory", "\n")

colnames(geno.bim) <- c("chr", "snp", "cM", "Position", "minor_allele", "major_allele")

setnames(cis.eqtl, "snps", "snp")
setnames(cis.eqtl, "SNPpos", "position")

cis.eqtl <- merge(cis.eqtl, geno.bim, by="snp")

cis.eqtl[, chr := replace(chr.y, chr.y == 23, "X")]
cis.eqtl[, varbeta := se^2]

cis.eqtl.loci <- split(cis.eqtl[,c("snp", "position", "beta", "varbeta", "chr", "minor_allele", "major_allele")], cis.eqtl$gene)

cat("Cis-eQTL Data Separated by Locus", "\n")

for (locus in names(cis.eqtl.loci)[1:100]) {
  fwrite(cis.eqtl.loci[[locus]], paste0(locus, ".csv"), quote=FALSE)
}

cat("Finished Writing", "\n")
