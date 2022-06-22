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

#----------------------------------------------------------
# Load Data
#----------------------------------------------------------

# Summary statistics from initial pass in cis regions
cis.pqtl <- readRDS("/nfs/users/nfs_n/nm18/gains_team282/proteomics/pqtl/pqtl_ms2019/cis_pqtl_all.RDS")

# Summary statistics from initial pass in trans regions
trans.pqtl <- readRDS("/nfs/users/nfs_n/nm18/gains_team282/proteomics/pqtl/pqtl_ms2019/trans_pqtl_all.RDS")

# SNP information
geno.bim <- fread("/nfs/users/nfs_n/nm18/gains_team282/Genotyping/All_genotyping_merged_filtered_b38_refiltered_rsID.bim", sep="\t")
colnames(geno.bim) <- c("chr", "snp", "cM", "Position", "minor_allele", "major_allele")

# Merge summary statistics with SNP information in cis regions
cis.pqtl <- merge(cis.pqtl, geno.bim, by.x="SNP", by.y="snp")
cis.pqtl$chr <- NULL
cis.pqtl$Position.y <- NULL
setnames(cis.pqtl, "Position.x", "Position")

# Merge summary statistics with SNP information in trans region
trans.pqtl <- merge(trans.pqtl, geno.bim, by.x="SNP", by.y="snp")
trans.pqtl$chr <- NULL
trans.pqtl$Position.y <- NULL
setnames(trans.pqtl, "Position.x", "Position")

# Split summary statistics by locus
cis.pqtl.loci <- split(
    cis.pqtl[,c("SNP", "Position", "pQTL_beta", "pQTL_SE", "Chr", "minor_allele", "major_allele")],
    paste0(cis.pqtl$Gene, "-", cis.pqtl$pQTL_Protein)
)
trans.pqtl.loci <- split(
    trans.pqtl[,c("SNP", "Position", "pQTL_beta", "pQTL_SE", "Chr", "minor_allele", "major_allele")],
    paste0(trans.pqtl$pQTL_Protein, "-", trans.pqtl$pQTL_Locus)
)

# Load genotypes for the relevant loci
cis.geno <- fread("/nfs/users/nfs_n/nm18/gains_team282/nikhil/data/genotypes/cis_pqtl_genotypes.raw", sep=" ", drop=2:6)
cis.geno <- as.data.frame(cis.geno)
colnames(cis.geno) <- gsub("X", "", colnames(cis.geno))
colnames(cis.geno) <- substr(colnames(cis.geno), 1, nchar(colnames(cis.geno)) - 2)
rownames(cis.geno) <- gsub("^GA", "", cis.geno[, 1])
cis.geno[, 1] <- NULL

trans.geno <- fread("/nfs/users/nfs_n/nm18/gains_team282/nikhil/data/genotypes/trans_pqtl_genotypes.raw", sep=" ", drop=2:6)
trans.geno <- as.data.frame(trans.geno)
colnames(trans.geno) <- gsub("X", "", colnames(trans.geno))
colnames(trans.geno) <- substr(colnames(trans.geno), 1, nchar(colnames(trans.geno)) - 2)
rownames(trans.geno) <- gsub("^GA", "", trans.geno[, 1])
trans.geno[, 1] <- NULL

# Load protein expression
protein.exp <- read.csv("/nfs/users/nfs_n/nm18/gains_team282/proteomics/MS2019_processed_data/data_291x1860_MS2019.csv", row.names=1)

protein.sample.info <- read.csv("/nfs/users/nfs_n/nm18/gains_team282/proteomics/MS2019_processed_data/sample_info_1860_MS2019.csv") %>%
    dplyr::mutate(Patient = gsub("^GA", "", Patient)) %>%
    dplyr::filter(Patient %in% rownames(cis.geno))

protein.exp <- protein.exp[, paste0("X", protein.sample.info$Injection, ".Razor.Intensity")]

n.samples <- nrow(protein.sample.info)
saveRDS(n.samples, "n_samples.RDS")

#----------------------------------------------------------
# Split Loci in Chromosome
#----------------------------------------------------------

# Only perform fine mapping for significant cis-eQTL
doParallel::registerDoParallel(cores=16)
foreach(locus=names(cis.pqtl.loci)) %dopar% {

    gene = gsub("-.*$", "", locus)
    protein = gsub("^.*-", "", locus)

    pqtl.locus <- cis.pqtl.loci[[locus]] %>%
        dplyr::arrange(Position) %>%
        dplyr::select(SNP, Beta=pQTL_beta, SE=pQTL_SE)

    pqtl.locus.geno <- cis.geno[, pqtl.locus$SNP]

    pqtl.locus.var.y <- var(as.numeric(protein.exp[protein, ]), na.rm=T)

    saveRDS(pqtl.locus, paste0(gene, ".summary.RDS"))
    saveRDS(pqtl.locus.geno, paste0(gene, ".genotypes.RDS"))
    saveRDS(pqtl.locus.var.y, paste0(gene, ".var.y.RDS"))
}

# Only perform fine mapping for significant cis-eQTL
doParallel::registerDoParallel(cores=16)
foreach(locus=names(trans.pqtl.loci)) %dopar% {

    protein = gsub("-.*$", "", locus)

    pqtl.locus <- trans.pqtl.loci[[locus]] %>%
        dplyr::arrange(Position) %>%
        dplyr::select(SNP, Beta=pQTL_beta, SE=pQTL_SE)

    pqtl.locus.geno <- trans.geno[, pqtl.locus$SNP]

    pqtl.locus.var.y <- var(as.numeric(protein.exp[protein, ]), na.rm=T)

    saveRDS(pqtl.locus, paste0(locus, ".summary.RDS"))
    saveRDS(pqtl.locus.geno, paste0(locus, ".genotypes.RDS"))
    saveRDS(pqtl.locus.var.y, paste0(locus, ".var.y.RDS"))
}
