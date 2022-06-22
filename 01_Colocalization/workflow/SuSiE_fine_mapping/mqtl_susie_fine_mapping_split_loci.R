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

module.ss.dir = "/nfs/users/nfs_n/nm18/gains_team282/nikhil/expression/eigengene_sva/wgcna_summary_statistics/"
module.qtl <- do.call(rbind, lapply(list.files(module.ss.dir, pattern="ME_[0-9]+_[0-9]+-.*\\.tsv"), function(file.name) {

    fread(paste0(module.ss.dir, file.name)) %>%
    as.data.frame() %>%
    dplyr::select(snp=1, beta=2, se=3, t=4, p=5) %>%
    dplyr::mutate(module.qtl=gsub("\\.tsv", "", file.name)) %>%
    dplyr::mutate(module=gsub("_[0-9]+-.*$", "", module.qtl)) %>%
    dplyr::mutate(pc=gsub("-.*$", "", gsub("ME_[0-9]+_", "", module.qtl))) %>%
    dplyr::mutate(qtl.locus=gsub("ME_[0-9]+_[0-9]+-", "", module.qtl)) %>%
    dplyr::mutate(qtl.locus.chr=gsub("\\:.*", "", qtl.locus)) %>%
    dplyr::mutate(qtl.locus.start=as.numeric(gsub(".*\\:", "", gsub("-.*$", "", qtl.locus)))) %>%
    dplyr::mutate(qtl.locus.end=as.numeric(gsub(".*-", "", qtl.locus)))
}))

# SNP information
geno.bim <- fread("/nfs/users/nfs_n/nm18/gains_team282/Genotyping/All_genotyping_merged_filtered_b38_refiltered_rsID.bim", sep="\t")
colnames(geno.bim) <- c("chr", "snp", "cM", "Position", "minor_allele", "major_allele")

# Merge summary statistics with SNP information
module.qtl <- merge(module.qtl, geno.bim, by="snp")

# Split summary statistics by locus
mqtl.loci <- split(
    module.qtl[,c("snp", "Position", "beta", "se", "chr", "minor_allele", "major_allele")],
    module.qtl$module.qtl
)

# Load genotypes for the relevant loci
geno <- fread("/nfs/users/nfs_n/nm18/gains_team282/nikhil/data/genotypes/eigengene_sva_ss_genotypes.raw", sep=" ", drop=2:6)
geno <- as.data.frame(geno)
colnames(geno) <- gsub("X", "", colnames(geno))
colnames(geno) <- substr(colnames(geno), 1, nchar(colnames(geno)) - 2)
rownames(geno) <- gsub("^GA", "", geno[, 1])
geno[, 1] <- NULL

# Load eigengene expression
eigengenes <- read.csv("/nfs/users/nfs_n/nm18/gains_team282/nikhil/expression/gene_expression/eigengenes.multiple.csv", row.names=1)

# Save number of samples
n.samples = nrow(eigengenes)
saveRDS(n.samples, "n_samples.RDS")

#----------------------------------------------------------
# Split Loci in Chromosome
#----------------------------------------------------------

# Only perform fine mapping for significant cis-eQTL
doParallel::registerDoParallel(cores=16)
foreach(locus=names(mqtl.loci)) %dopar% {

    me = gsub("-.*$", "", locus)

    mqtl.locus <- mqtl.loci[[locus]] %>%
        dplyr::arrange(Position) %>%
        dplyr::select(SNP=snp, Beta=beta, SE=se)

    mqtl.locus.geno <- geno[, mqtl.locus$SNP]

    mqtl.locus.var.y <- var(eigengenes[, me], na.rm=T)

    saveRDS(mqtl.locus, paste0(locus, ".summary.RDS"))
    saveRDS(mqtl.locus.geno, paste0(locus, ".genotypes.RDS"))
    saveRDS(mqtl.locus.var.y, paste0(locus, ".var.y.RDS"))
}
