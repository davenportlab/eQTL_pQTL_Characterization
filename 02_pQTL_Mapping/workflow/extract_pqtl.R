#----------------------------------------------------------
# Extract pQTL from Summary Results
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
library(GenomicRanges)

options(stringsAsFactors = FALSE)

args = commandArgs(trailingOnly = TRUE)
protein.input = args[1]
n.proteins = as.numeric(args[2])

#----------------------------------------------------------
# Load Data
#----------------------------------------------------------

# eQTL-pQTL Colocalization Metadata
metadata <- read.table("/nfs/users/nfs_n/nm18/gains_team282/nikhil/colocalization/eQTL_pQTL_metadata.tsv", sep="\t", header=T)

# Summary statistics from pQTL mapping
summary <- do.call(rbind, mclapply(1:4, function(i) {
    readRDS(paste0("summary_stats/", protein.input, "_group_", i, ".rds"))
}, mc.cores=4)) %>%
    as.data.frame() %>%
    dplyr::mutate(SNP=substr(SNP, 1, nchar(SNP) - 2)) %>%
    dplyr::mutate(pQTL_beta=as.numeric(pQTL_beta)) %>%
    dplyr::mutate(pQTL_SE=as.numeric(pQTL_SE)) %>%
    dplyr::mutate(pQTL_t=as.numeric(pQTL_t)) %>%
    dplyr::mutate(pQTL_pval=as.numeric(pQTL_pval))
    
# SNP information
geno.bim <- fread("/nfs/users/nfs_n/nm18/gains_team282/Genotyping/All_genotyping_merged_filtered_b38_refiltered_rsID.bim", sep="\t") %>%
    as.data.frame()
colnames(geno.bim) <- c("Chr", "SNP", "cM", "Position", "minor_allele", "major_allele")

# Integrate SNP information with summary statistics
summary.full <- merge(summary, geno.bim, by="SNP") %>%
    dplyr::select(SNP, Chr, Position, pQTL_beta, pQTL_SE, pQTL_t, pQTL_pval, pQTL_Protein=Protein) %>%
    dplyr::arrange(Chr, Position)

#----------------------------------------------------------
# Identify pQTL
#----------------------------------------------------------

threshold <- 0.05 / (n.proteins * nrow(summary.full))
window.size <- 1e6

pqtl.loci <- summary.full[summary.full$pQTL_pval < threshold,] %>%
    dplyr::mutate(start=Position - (window.size / 2) + 1) %>%
    dplyr::mutate(end=Position + (window.size / 2))

if (nrow(pqtl.loci) > 0) {

    pqtl.loci <- pqtl.loci %>%
        makeGRangesFromDataFrame(keep.extra.columns=TRUE) %>%
        reduce() %>%
        as.data.frame()

    pqtl.summary <- do.call(rbind, mclapply(1:nrow(pqtl.loci), function(i) {

        pqtl.locus = paste0("pQTL_", i)

        summary.full %>%
            dplyr::filter(Chr == pqtl.loci$seqnames[i]) %>%
            dplyr::filter(Position >= pqtl.loci$start[i]) %>%
            dplyr::filter(Position <= pqtl.loci$end[i]) %>%
            dplyr::mutate(pQTL_Locus=pqtl.locus)
    }, mc.cores=16))
} else {

    pqtl.summary <- data.frame()
}

saveRDS(pqtl.summary, "pqtl.summary.RDS")

#----------------------------------------------------------
# Extract cis-pQTL Region
#----------------------------------------------------------

protein.metadata <- metadata %>%
    dplyr::filter(Accession == summary$Protein[1]) %>%
    dplyr::filter(!is.na(tss))

if (nrow(protein.metadata) > 0) {
    
    window.size <- 1e6

    cis.pqtl.summary <- do.call(rbind, mclapply(1:nrow(protein.metadata), function(i) {

        summary.full %>%
            dplyr::filter(Chr == protein.metadata$seqname[i]) %>%
            dplyr::filter(Position >= protein.metadata$tss[i] - (window.size / 2) + 1) %>%
            dplyr::filter(Position <= protein.metadata$tss[i] + (window.size / 2)) %>%
            dplyr::mutate(Gene = protein.metadata$gene_id[i])
    }, mc.cores=16))
} else {

    cis.pqtl.summary <- data.frame()
}

saveRDS(cis.pqtl.summary, "cis.pqtl.summary.RDS")
