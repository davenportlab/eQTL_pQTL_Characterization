#----------------------------------------------------------
# Prepare Loci for Conditional Effects from FINEMAP
# Created: 20 March 2022
#----------------------------------------------------------

args = commandArgs(trailingOnly=TRUE)
chr = args[1]

#----------------------------------------------------------
# Environment
#----------------------------------------------------------

library(tidyverse)
library(GenomicRanges)
library(data.table)
library(parallel)
library(foreach)

# SNP information
geno.bim <- fread("/nfs/users/nfs_n/nm18/gains_team282/Genotyping/All_genotyping_merged_filtered_b38_refiltered_rsID.bim", sep="\t")
colnames(geno.bim) <- c("chr", "snp", "cM", "Position", "minor_allele", "major_allele")

# Add allele frequencies
allele.freqs <- fread("plink.frq", header=TRUE)
geno.bim <- merge(geno.bim, allele.freqs[,c("SNP", "MAF", "NCHROBS")], by.x="snp", by.y="SNP")

# FINEMAP Analysis Results
cs <- fread(paste0("/nfs/users/nfs_n/nm18/gains_team282/nikhil/colocalization/cis_eqtl/fine_mapping/FINEMAP/chr", chr, "_credible_sets.tsv")) %>%
    as.data.frame() %>%
    dplyr::group_by(Gene) %>%
    dplyr::slice_max(Post_Prob_k, n=1, with_ties=TRUE) %>%
    dplyr::group_by(Gene, k, Credibility_Set) %>%
    dplyr::slice_max(Prob_SNP_in_CS, n=1, with_ties=FALSE)

cs.loci <- split(cs, cs$Gene)

# Summary statistics from initial pass
cis.eqtl.summary <- readRDS("/nfs/users/nfs_n/nm18/gains_team282/eqtl/cisresults/ciseqtl_all.rds")
cis.eqtl.summary <- cis.eqtl.summary[cis.eqtl.summary$gene %in% cs$Gene,]
cis.eqtl.summary <- merge(cis.eqtl.summary, geno.bim, by.x="snps", by.y="snp")
cis.eqtl.summary$chr.y <- NULL
setnames(cis.eqtl.summary, "chr.x", "chr")
cis.eqtl.loci <- split(cis.eqtl.summary, cis.eqtl.summary$gene)

# Number of signals at each locus
cs.n.signals <- cs %>%
    dplyr::group_by(Gene) %>%
    dplyr::summarize(Signals=n())

#----------------------------------------------------------
# Single Signal Summary Statistics
#----------------------------------------------------------

single.signal.stats <- do.call(rbind, mclapply(cs.n.signals$Gene[cs.n.signals$Signals == 1], function(locus) {

    cis.eqtl.loci[[locus]] %>%
        dplyr::mutate(Gene=locus, Signal=1) %>%
        dplyr::select(Gene, Signal, chr, snps, SNPpos, beta, se, pvalue)
}, mc.cores=16))

write.table(single.signal.stats, "single_signal_stats.txt", sep="\t", row.names=F, col.names=F, quote=F)

#----------------------------------------------------------
# Multiple Signal Summary Statistics
#----------------------------------------------------------

# Only identify conditional summary statistics for loci with more than one signal

doParallel::registerDoParallel(cores=16)
foreach(locus=cs.n.signals$Gene[cs.n.signals$Signals > 1]) %dopar% {

    eqtl.locus <- cis.eqtl.loci[[locus]] %>%
        dplyr::mutate(N = NCHROBS / 2) %>%
        dplyr::select(SNP=snps, A1=minor_allele, A2=major_allele, freq=MAF, b=beta, se, p=pvalue, N)

    write.table(eqtl.locus, paste0(locus, ".ma"), sep=" ", quote=F, row.names=F)

    # Save SNPs at eQTL locus
    write.table(eqtl.locus$SNP, paste0(locus, ".snps"), row.names=F, col.names=F, quote=F)

    # Save conditional SNP list for each conditional signal
    conditional.snps <- cs.loci[[locus]]$SNP
    for (i in 1:length(conditional.snps)) {
        snp = conditional.snps[i]
        number = cs.loci[[locus]]$Credibility_Set[i]
        write.table(setdiff(conditional.snps, snp), paste0(locus, "-", number, ".cond.snps"), row.names=F, col.names=F, quote=F)
    }
}
