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

# Significant cis-eQTL used to filter the loci tested using FINEMAP
cis.eqtl <- readRDS("/nfs/users/nfs_n/nm18/gains_team282/eqtl/cisresults/cisqtl_all_significant.rds")

# Conditional summary statistics
cis.eqtl.summary.conditional <- fread(paste0("/nfs/users/nfs_n/nm18/gains_team282/nikhil/colocalization/cis_eqtl/conditional_effects/LD/chr", chr.input, "_conditional_cis_eQTL_summary_statistics.tsv"))

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

    master.file <- data.frame(
        z=paste0(locus, ".z"),
        ld=paste0(locus, ".ld"),
        snp=paste0(locus, ".snp"),
        config=paste0(locus, ".config"),
        cred=paste0(locus, ".cred"),
        log=paste0(locus, ".log"),
        n_samples=n.samples
    )

    write.table(z.file, paste0("full/", locus, ".z"), sep=" ", quote=F, row.names=F)
    write.table(ld.file, paste0("full/", locus, ".ld"), sep=" ", quote=F, row.names=F, col.names=F)
    write.table(master.file, paste0("full/", locus, ".master"), sep=";", quote=F, row.names=F)
}

# Only perform fine mapping for significant cis-eQTL
doParallel::registerDoParallel(cores=16)
foreach(locus=names(cis.eqtl.conditional.loci)) %dopar% {

    z.file <- cis.eqtl.conditional.loci[[locus]] %>%
        dplyr::mutate(maf=maf[SNP,1]) %>%
        dplyr::select(rsid=SNP, chromosome=Chr, position=Position, allele1=minor_allele, allele2=major_allele, maf, beta=Beta, se=SE) %>%
        dplyr::filter(!is.na(beta))

    eqtl.locus.geno <- chr.geno[, z.file$rsid]
    ld.file <- cor(eqtl.locus.geno, use="pairwise.complete.obs")

    master.file <- data.frame(
        z=paste0(locus, ".z"),
        ld=paste0(locus, ".ld"),
        snp=paste0(locus, ".snp"),
        config=paste0(locus, ".config"),
        cred=paste0(locus, ".cred"),
        log=paste0(locus, ".log"),
        n_samples=n.samples
    )

    write.table(z.file, paste0("conditional/", locus, ".z"), sep=" ", quote=F, row.names=F)
    write.table(ld.file, paste0("conditional/", locus, ".ld"), sep=" ", quote=F, row.names=F, col.names=F)
    write.table(master.file, paste0("conditional/", locus, ".master"), sep=";", quote=F, row.names=F)
}
