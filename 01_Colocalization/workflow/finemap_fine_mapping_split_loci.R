#----------------------------------------------------------
# Fine Mapping with FINEMAP
# Created: 17 February 2022
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

# Split summary statistics by locus
cis.eqtl.loci <- split(cis.eqtl.summary[,c("snp", "position", "beta", "se", "chr", "minor_allele", "major_allele")], cis.eqtl.summary$gene)

# Load genotypes for the relevant loci
chr.geno <- fread(paste0("/nfs/users/nfs_n/nm18/gains_team282/nikhil/data/genotypes/eqtl_genotypes_", chr.input, ".raw"), sep=" ", drop=2:6)
chr.geno <- as.data.frame(chr.geno)
colnames(chr.geno) <- gsub("X", "", colnames(chr.geno))
colnames(chr.geno) <- sapply(strsplit(colnames(chr.geno), "_"), function(x) { x[1] })
rownames(chr.geno) <- gsub("^GA", "", chr.geno[, 1])
chr.geno[, 1] <- NULL

# Calculate Minor Allele Frequencies (MAFs)
maf <- as.matrix(
    apply(
        chr.geno, 
        2, 
        function(x) { sum(x, na.rm=T) }
    ) / (2 * nrow(chr.geno))
)

# Load gene expression to get number of samples
gene.exp <- read.table("/lustre/scratch119/humgen/projects/gains_team282/eqtl/data/logcpm_864_20412_hla.txt")
n.samples <- ncol(gene.exp)

# Save master file for FINEMAP
master.file <- data.frame(Locus=names(cis.eqtl.loci)) %>% 
    dplyr::mutate(z=paste0(Locus, ".z")) %>%
    dplyr::mutate(ld=paste0(Locus, ".ld")) %>%
    dplyr::mutate(snp=paste0(Locus, ".snp")) %>%
    dplyr::mutate(config=paste0(Locus, ".config")) %>%
    dplyr::mutate(cred=paste0(Locus, ".cred")) %>%
    dplyr::mutate(log=paste0(Locus, ".log")) %>%
    dplyr::mutate(n_samples=n.samples) %>%
    dplyr::select(z, ld, snp, config, cred, log, n_samples)

write.table(master.file, "master", sep=";", quote=F, row.names=F)

#----------------------------------------------------------
# Split Loci in Chromosome
#----------------------------------------------------------

# Only perform fine mapping for significant cis-eQTL
doParallel::registerDoParallel(cores=16)
foreach(locus=names(cis.eqtl.loci)) %dopar% {

    z.file <- cis.eqtl.loci[[locus]] %>%
        dplyr::mutate(maf=maf[snp,1]) %>%
        dplyr::select(rsid=snp, chromosome=chr, position, allele1=minor_allele, allele2=major_allele, maf, beta, se)

    eqtl.locus.geno <- chr.geno[, z.file$rsid]
    ld.file <- cor(eqtl.locus.geno, use="pairwise.complete.obs")

    write.table(z.file, paste0(locus, ".z"), sep=" ", quote=F, row.names=F)
    write.table(ld.file, paste0(locus, ".ld"), sep=" ", quote=F, row.names=F, col.names=F)
}
