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
library(parallel)
library(foreach)

options(stringsAsFactors = FALSE)

args = commandArgs(trailingOnly = TRUE)
chr.input = args[1]

#----------------------------------------------------------
# Load Data
#----------------------------------------------------------

# Summary statistics from initial pass
cis.eqtl.summary <- readRDS("/nfs/users/nfs_n/nm18/gains_team282/eqtl/cisresults/ciseqtl_all.rds")

# Significant cis-eQTL used to filter the loci tested using SuSiE
cis.eqtl <- readRDS("/nfs/users/nfs_n/nm18/gains_team282/eqtl/cisresults/cisqtl_all_significant.rds")

# Conditional summary statistics
cis.eqtl.summary.conditional <- fread(paste0("/nfs/users/nfs_n/nm18/gains_team282/nikhil/colocalization/cis_eqtl/conditional_effects/LMM/chr", chr.input, "_conditional_cis_eQTL_summary_statistics.tsv"))

# SNP information
geno.bim <- fread("/nfs/users/nfs_n/nm18/gains_team282/Genotyping/All_genotyping_merged_filtered_b38_refiltered_rsID.bim", sep="\t")
colnames(geno.bim) <- c("chr", "snp", "cM", "Position", "minor_allele", "major_allele")

# Rename columns in the summary statistics
setnames(cis.eqtl.summary, "snps", "snp")
setnames(cis.eqtl.summary, "SNPpos", "position")

# Subset summary statistics based on loci with significant cis-eQTL
cis.eqtl.summary <- subset(cis.eqtl.summary, chr == chr.input & gene %in% cis.eqtl$gene)

# Merge summary statistics with SNP information
cis.eqtl.summary <- merge(cis.eqtl.summary, geno.bim, by="snp")
cis.eqtl.summary$chr.y <- NULL
setnames(cis.eqtl.summary, "chr.x", "chr")

cis.eqtl.summary.conditional <- merge(cis.eqtl.summary.conditional, geno.bim, by.x="SNP", by.y="snp")
cis.eqtl.summary.conditional$Position.y <- NULL
setnames(cis.eqtl.summary.conditional, "Position.x", "Position")

# Filter cis-eQTL to those that are significant
cis.eqtl.summary <- cis.eqtl.summary %>%
    dplyr::filter(gene %in% cis.eqtl$gene)

# Split summary statistics by locus
cis.eqtl.loci <- split(cis.eqtl.summary[,c("snp", "position", "beta", "se", "chr", "minor_allele", "major_allele")], cis.eqtl.summary$gene)

cis.eqtl.conditional.loci <- split(cis.eqtl.summary.conditional, paste0(cis.eqtl.summary.conditional$Gene, "-", cis.eqtl.summary.conditional$Signal))

# Load genotypes for the relevant loci
chr.geno <- fread(paste0("/nfs/users/nfs_n/nm18/gains_team282/nikhil/data/genotypes/eqtl_genotypes_", chr.input, ".raw"), sep=" ", drop=2:6)
chr.geno <- as.data.frame(chr.geno)
colnames(chr.geno) <- gsub("X", "", colnames(chr.geno))
colnames(chr.geno) <- sapply(strsplit(colnames(chr.geno), "_"), function(x) { x[1] })
rownames(chr.geno) <- gsub("^GA", "", chr.geno[, 1])
chr.geno[, 1] <- NULL

# Load gene expression
gene.exp <- read.table("/lustre/scratch119/humgen/projects/gains_team282/eqtl/data/logcpm_864_20412_hla.txt")
gene.exp <- t(gene.exp)
rownames(gene.exp) <- gsub("^GA", "", rownames(gene.exp))

n.samples <- length(unique(sapply(strsplit(rownames(gene.exp), "_"), function(x) { x[1] })))
saveRDS(n.samples, "n_samples.RDS")

#----------------------------------------------------------
# Split Loci in Chromosome
#----------------------------------------------------------

# Only perform fine mapping for significant cis-eQTL
doParallel::registerDoParallel(cores=16)
foreach(locus=names(cis.eqtl.loci)) %dopar% {

    eqtl.locus <- cis.eqtl.loci[[locus]] %>%
        dplyr::arrange(position)

    eqtl.locus.geno <- chr.geno[, eqtl.locus$snp]

    eqtl.locus.var.y <- var(as.numeric(gene.exp[, locus]), na.rm=T)

    saveRDS(eqtl.locus, paste0("full/", locus, ".summary.RDS"))
    saveRDS(eqtl.locus.geno, paste0("full/", locus, ".genotypes.RDS"))
    saveRDS(eqtl.locus.var.y, paste0("full/", locus, ".var.y.RDS"))
}

# Only perform fine mapping for significant cis-eQTL
doParallel::registerDoParallel(cores=16)
foreach(locus=names(cis.eqtl.conditional.loci)) %dopar% {

    eqtl.locus <- cis.eqtl.conditional.loci[[locus]] %>%
        dplyr::arrange(Position) %>%
        dplyr::select(snp=SNP, beta=Beta, se=SE)

    eqtl.locus.geno <- chr.geno[, eqtl.locus$snp]

    gene.name = gsub("-.*$", "", locus)
    eqtl.locus.var.y <- var(as.numeric(gene.exp[, gene.name]), na.rm=T)

    saveRDS(eqtl.locus, paste0("conditional/", locus, ".summary.RDS"))
    saveRDS(eqtl.locus.geno, paste0("conditional/", locus, ".genotypes.RDS"))
    saveRDS(eqtl.locus.var.y, paste0("conditional/", locus, ".var.y.RDS"))
}
