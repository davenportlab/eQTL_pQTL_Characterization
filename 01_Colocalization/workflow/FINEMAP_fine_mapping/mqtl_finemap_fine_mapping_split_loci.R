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

# Calculate Minor Allele Frequencies (MAFs)
maf <- as.matrix(
    apply(
        geno, 
        2, 
        function(x) { sum(x, na.rm=T) }
    ) / (2 * nrow(geno))
)

# Load gene expression to get number of samples
n.samples = system("head -n 1 /lustre/scratch119/humgen/projects/gains_team282/eqtl/data/logcpm_864_20412_hla.txt | sed 's/\\t/\\n/g' | wc -l", intern=TRUE)
n.samples = as.numeric(n.samples)

#----------------------------------------------------------
# Write Out Module QTL Information
#----------------------------------------------------------

doParallel::registerDoParallel(cores=16)
foreach(locus=names(mqtl.loci)) %dopar% {

    z.file <- mqtl.loci[[locus]] %>%
        dplyr::mutate(maf=maf[snp,1]) %>%
        dplyr::select(rsid=snp, chromosome=chr, position=Position, allele1=minor_allele, allele2=major_allele, maf, beta=beta, se=se)

    locus.geno <- geno[, z.file$rsid]
    ld.file <- cor(locus.geno, use="pairwise.complete.obs")

    master.file <- data.frame(
        z=paste0(locus, ".z"),
        ld=paste0(locus, ".ld"),
        snp=paste0(locus, ".snp"),
        config=paste0(locus, ".config"),
        cred=paste0(locus, ".cred"),
        log=paste0(locus, ".log"),
        n_samples=n.samples
    )

    write.table(z.file, paste0(locus, ".z"), sep=" ", quote=F, row.names=F)
    write.table(ld.file, paste0(locus, ".ld"), sep=" ", quote=F, row.names=F, col.names=F)
    write.table(master.file, paste0(locus, ".master"), sep=";", quote=F, row.names=F)
}
