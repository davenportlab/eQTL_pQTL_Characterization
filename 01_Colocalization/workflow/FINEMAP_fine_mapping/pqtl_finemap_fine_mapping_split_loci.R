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
    cis.pqtl$Gene
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

# Calculate Minor Allele Frequencies (MAFs)
cis.maf <- as.matrix(
    apply(
        cis.geno, 
        2, 
        function(x) { sum(x, na.rm=T) }
    ) / (2 * nrow(cis.geno))
)

trans.maf <- as.matrix(
    apply(
        trans.geno,
        2,
        function(x) { sum(x, na.rm=T) }
    ) / (2 * nrow(trans.geno))
)

# Load protein expression to get number of samples
protein.sample.info <- read.csv("/nfs/users/nfs_n/nm18/gains_team282/proteomics/MS2019_processed_data/sample_info_1860_MS2019.csv") %>%
    dplyr::mutate(Patient = gsub("^GA", "", Patient)) %>%
    dplyr::filter(Patient %in% rownames(cis.geno))
n.samples <- nrow(protein.sample.info)

#----------------------------------------------------------
# Write Out Cis-pQTL Information
#----------------------------------------------------------

dir.create("cis_pqtl/")

doParallel::registerDoParallel(cores=16)
foreach(locus=names(cis.pqtl.loci)) %dopar% {

    z.file <- cis.pqtl.loci[[locus]] %>%
        dplyr::mutate(maf=cis.maf[SNP,1]) %>%
        dplyr::select(rsid=SNP, chromosome=Chr, position=Position, allele1=minor_allele, allele2=major_allele, maf, beta=pQTL_beta, se=pQTL_SE)

    pqtl.locus.geno <- cis.geno[, z.file$rsid]
    ld.file <- cor(pqtl.locus.geno, use="pairwise.complete.obs")

    master.file <- data.frame(
        z=paste0(locus, ".z"),
        ld=paste0(locus, ".ld"),
        snp=paste0(locus, ".snp"),
        config=paste0(locus, ".config"),
        cred=paste0(locus, ".cred"),
        log=paste0(locus, ".log"),
        n_samples=n.samples
    )

    write.table(z.file, paste0("cis_pqtl/", locus, ".z"), sep=" ", quote=F, row.names=F)
    write.table(ld.file, paste0("cis_pqtl/", locus, ".ld"), sep=" ", quote=F, row.names=F, col.names=F)
    write.table(master.file, paste0("cis_pqtl/", locus, ".master"), sep=";", quote=F, row.names=F)
}

#----------------------------------------------------------
# Write Out Trans-pQTL Information
#----------------------------------------------------------

dir.create("trans_pqtl/")

doParallel::registerDoParallel(cores=16)
foreach(locus=names(trans.pqtl.loci)) %dopar% {

    z.file <- trans.pqtl.loci[[locus]] %>%
        dplyr::mutate(maf=trans.maf[SNP,1]) %>%
        dplyr::select(rsid=SNP, chromosome=Chr, position=Position, allele1=minor_allele, allele2=major_allele, maf, beta=pQTL_beta, se=pQTL_SE)

    pqtl.locus.geno <- trans.geno[, z.file$rsid]
    ld.file <- cor(pqtl.locus.geno, use="pairwise.complete.obs")

    master.file <- data.frame(
        z=paste0(locus, ".z"),
        ld=paste0(locus, ".ld"),
        snp=paste0(locus, ".snp"),
        config=paste0(locus, ".config"),
        cred=paste0(locus, ".cred"),
        log=paste0(locus, ".log"),
        n_samples=n.samples
    )

    write.table(z.file, paste0("trans_pqtl/", locus, ".z"), sep=" ", quote=F, row.names=F)
    write.table(ld.file, paste0("trans_pqtl/", locus, ".ld"), sep=" ", quote=F, row.names=F, col.names=F)
    write.table(master.file, paste0("trans_pqtl/", locus, ".master"), sep=";", quote=F, row.names=F)
}
