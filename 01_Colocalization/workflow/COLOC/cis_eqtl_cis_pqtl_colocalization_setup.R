#----------------------------------------------------------
# Cis-eQTL Cis-pQTL Overlap
# Created: 17 March 2022
#----------------------------------------------------------

#----------------------------------------------------------
# Environment
#----------------------------------------------------------

library(tidyverse)
library(GenomicRanges)
library(data.table)

#----------------------------------------------------------
# Read Cis-eQTL and Cis-pQTL
#----------------------------------------------------------

cis.eqtl <- read.table("/nfs/users/nfs_n/nm18/gains_team282/eqtl/cisresults/eigenMT/ciseqtl_eigenMT_corrected.txt") %>%
    dplyr::filter(Sig)

prot.info <- read.table("/nfs/users/nfs_n/nm18/gains_team282/nikhil/colocalization/eQTL_pQTL_metadata.tsv", sep="\t", header=TRUE)

common.genes <- intersect(cis.eqtl$gene, prot.info$Gene_ID)

common.genes.chrs <- cis.eqtl$chr[cis.eqtl$gene %in% common.genes]

common.genes.split <- split(common.genes, common.genes.chrs)

#----------------------------------------------------------
# Load Summary Information from Mapping
#----------------------------------------------------------

# SNP information
geno.bim <- fread("/nfs/users/nfs_n/nm18/gains_team282/Genotyping/All_genotyping_merged_filtered_b38_refiltered_rsID.bim", sep="\t")
colnames(geno.bim) <- c("chr", "snp", "cM", "Position", "minor_allele", "major_allele")

# Summary statistics from initial pass
cis.eqtl.summary <- readRDS("/nfs/users/nfs_n/nm18/gains_team282/eqtl/cisresults/ciseqtl_all.rds")
cis.pqtl.summary <- readRDS("/nfs/users/nfs_n/nm18/gains_team282/proteomics/pqtl/pqtl_ms2019/cis_pqtl_all.RDS")

cis.eqtl.summary <- subset(cis.eqtl.summary, cis.eqtl.summary$gene %in% common.genes)
cis.pqtl.summary <- subset(cis.pqtl.summary, cis.pqtl.summary$Gene %in% common.genes)

cis.eqtl.summary <- merge(cis.eqtl.summary, geno.bim, by.x="snps", by.y="snp")
cis.eqtl.summary$chr.y <- NULL
setnames(cis.eqtl.summary, "chr.x", "chr")

cis.pqtl.summary <- merge(cis.pqtl.summary, geno.bim, by.x="SNP", by.y="snp")
cis.pqtl.summary$Position.y <- NULL
setnames(cis.pqtl.summary, "Position.x", "Position")

# Split summary statistics by locus
cis.eqtl.loci <- split(cis.eqtl.summary[,c("snps", "Position", "beta", "se", "chr", "minor_allele", "major_allele")], cis.eqtl.summary$gene)
cis.pqtl.loci <- split(cis.pqtl.summary[,c("SNP", "Position", "pQTL_beta", "pQTL_SE", "chr", "minor_allele", "major_allele", "pQTL_Protein")], cis.pqtl.summary$Gene)

#----------------------------------------------------------
# Load Protein Genotypes
#----------------------------------------------------------

# Load pQTL genotypes
pqtl.geno <- fread("/nfs/users/nfs_n/nm18/gains_team282/nikhil/data/genotypes/cis_pqtl_genotypes.raw", sep=" ", drop=2:6)
pqtl.geno <- as.data.frame(pqtl.geno)
colnames(pqtl.geno) <- gsub("X", "", colnames(pqtl.geno))
colnames(pqtl.geno) <- sapply(strsplit(colnames(pqtl.geno), "_"), function(x) { x[1] })
rownames(pqtl.geno) <- gsub("^GA", "", pqtl.geno[, 1])
pqtl.geno[, 1] <- NULL

#----------------------------------------------------------
# Load Gene and Protein Expression
#----------------------------------------------------------

gene.exp <- read.table("/lustre/scratch119/humgen/projects/gains_team282/eqtl/data/logcpm_864_20412_hla.txt")
gene.exp <- t(gene.exp)
rownames(gene.exp) <- gsub("^GA", "", rownames(gene.exp))

protein.exp <- read.csv("/lustre/scratch119/humgen/projects/gains_team282/proteomics/MS2019_processed_data/data_291x1860_MS2019.csv", row.names=1)
protein.exp <- t(protein.exp)
rownames(protein.exp) <- gsub("^X", "", gsub("\\..*$", "", rownames(protein.exp)))
protein.info <- read.csv("/lustre/scratch119/humgen/projects/gains_team282/proteomics/MS2019_processed_data/sample_info_1860_MS2019.csv", row.names=1)
protein.info$Patient <- gsub("^GA", "", protein.info$Patient)
protein.info <- protein.info %>%
    dplyr::filter(Patient %in% rownames(pqtl.geno))
protein.exp <- protein.exp[paste0(protein.info$Injection),]

result <- do.call(rbind, lapply(unique(common.genes.chrs), function(chr) {

    # Load genotypes for the relevant loci
    chr.geno <- fread(paste0("/nfs/users/nfs_n/nm18/gains_team282/nikhil/data/genotypes/eqtl_genotypes_", chr, ".raw"), sep=" ", drop=2:6)
    chr.geno <- as.data.frame(chr.geno)
    colnames(chr.geno) <- gsub("X", "", colnames(chr.geno))
    colnames(chr.geno) <- sapply(strsplit(colnames(chr.geno), "_"), function(x) { x[1] })
    rownames(chr.geno) <- gsub("^GA", "", chr.geno[, 1])
    chr.geno[, 1] <- NULL

    do.call(rbind, lapply(common.genes.split[[paste0(chr)]], function(gene) {

        cis.eqtl.locus <- list()
        cis.eqtl.locus$beta <- cis.eqtl.loci[[gene]]$beta
        cis.eqtl.locus$varbeta <- cis.eqtl.loci[[gene]]$se^2
        cis.eqtl.locus$snp <- cis.eqtl.loci[[gene]]$snps
        cis.eqtl.locus$position <- cis.eqtl.loci[[gene]]$Position
        cis.eqtl.locus$type <- "quant"
        cis.eqtl.locus$sdY <- sd(gene.exp[,gene], na.rm=TRUE)
        cis.eqtl.locus$LD <- cor(chr.geno[,cis.eqtl.locus$snp], use="pairwise.complete.obs")

        cis.pqtl.locus <- list()
        cis.pqtl.locus$beta <- cis.pqtl.loci[[gene]]$pQTL_beta
        cis.pqtl.locus$varbeta <- cis.pqtl.loci[[gene]]$pQTL_SE^2
        cis.pqtl.locus$snp <- cis.pqtl.loci[[gene]]$SNP
        cis.pqtl.locus$position <- cis.pqtl.loci[[gene]]$Position
        cis.pqtl.locus$type <- "quant"
        cis.pqtl.locus$sdY <- sd(protein.exp[,cis.pqtl.loci[[gene]]$pQTL_Protein[1]], na.rm=TRUE)
        cis.pqtl.locus$LD <- cor(pqtl.geno[,cis.pqtl.locus$snp], use="pairwise.complete.obs")

        abf.res = coloc.abf(cis.eqtl.locus, cis.pqtl.locus)

        return(abf.res$summary)
    }))
}))
