#----------------------------------------------------------
# Prepare Loci for Conditional Effects from FINEMAP
# Created: 21 March 2022
#----------------------------------------------------------

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
cs.cis <- fread("/nfs/users/nfs_n/nm18/gains_team282/nikhil/colocalization/pqtl/fine_mapping/FINEMAP/cis_credible_sets.tsv") %>%
    as.data.frame() %>%
    dplyr::group_by(Gene) %>%
    dplyr::slice_max(Post_Prob_k, n=1, with_ties=TRUE) %>%
    dplyr::group_by(Gene, k, Credibility_Set) %>%
    dplyr::slice_max(Prob_SNP_in_CS, n=1, with_ties=FALSE)

cs.trans <- fread("/nfs/users/nfs_n/nm18/gains_team282/nikhil/colocalization/pqtl/fine_mapping/FINEMAP/trans_credible_sets.tsv") %>%
    as.data.frame() %>%
    dplyr::group_by(Gene) %>%
    dplyr::slice_max(Post_Prob_k, n=1, with_ties=TRUE) %>%
    dplyr::group_by(Gene, k, Credibility_Set) %>%
    dplyr::slice_max(Prob_SNP_in_CS, n=1, with_ties=FALSE)

cs.cis.loci <- split(cs.cis, cs.cis$Gene)
cs.trans.loci <- split(cs.trans, cs.trans$Gene)

# Summary statistics from initial pass
cis.pqtl.summary <- readRDS("/nfs/users/nfs_n/nm18/gains_team282/proteomics/pqtl/pqtl_ms2019/cis_pqtl_all.RDS")
cis.pqtl.summary <- cis.pqtl.summary[cis.pqtl.summary$Gene %in% cs.cis$Gene,]
cis.pqtl.summary <- merge(cis.pqtl.summary, geno.bim, by.x="SNP", by.y="snp")
cis.pqtl.summary$Position.y <- NULL
setnames(cis.pqtl.summary, "Position.x", "Position")
cis.pqtl.loci <- split(cis.pqtl.summary, cis.pqtl.summary$Gene)

trans.pqtl.summary <- readRDS("/nfs/users/nfs_n/nm18/gains_team282/proteomics/pqtl/pqtl_ms2019/trans_pqtl_all.RDS")
trans.pqtl.summary$ID <- with(trans.pqtl.summary, paste0(pQTL_Protein, "-", pQTL_Locus))
trans.pqtl.summary <- trans.pqtl.summary[trans.pqtl.summary$ID %in% cs.trans$Gene,]
trans.pqtl.summary <- merge(trans.pqtl.summary, geno.bim, by.x="SNP", by.y="snp")
trans.pqtl.summary$Position.y <- NULL
setnames(trans.pqtl.summary, "Position.x", "Position")
trans.pqtl.loci <- split(trans.pqtl.summary, trans.pqtl.summary$ID)

# Number of signals at each locus
cs.cis.n.signals <- cs.cis %>%
    dplyr::group_by(Gene) %>%
    dplyr::summarize(Signals=n())

cs.trans.n.signals <- cs.trans %>%
    dplyr::group_by(Gene) %>%
    dplyr::summarize(Signals=n())

#----------------------------------------------------------
# Single Signal Summary Statistics
#----------------------------------------------------------

# Cis-pQTL
cis.single.signal.stats <- do.call(rbind, mclapply(cs.cis.n.signals$Gene[cs.cis.n.signals$Signals == 1], function(locus) {

    cis.pqtl.loci[[locus]] %>%
        dplyr::mutate(Gene=locus, Signal=1) %>%
        dplyr::select(Gene, Signal, Chr, SNP, Position, pQTL_beta, pQTL_SE, pQTL_pval)
}, mc.cores=16))

write.table(cis.single.signal.stats, "cis_single_signal_stats.txt", sep="\t", row.names=F, col.names=F, quote=F)

# Trans-pQTL
trans.single.signal.stats <- do.call(rbind, mclapply(cs.trans.n.signals$Gene[cs.trans.n.signals$Signals == 1], function(locus) {

    trans.pqtl.loci[[locus]] %>%
        dplyr::mutate(Gene=locus, Signal=1) %>%
        dplyr::select(Gene, Signal, Chr, SNP, Position, pQTL_beta, pQTL_SE, pQTL_pval)
}, mc.cores=16))

write.table(trans.single.signal.stats, "trans_single_signal_stats.txt", sep="\t", row.names=F, col.names=F, quote=F)

#----------------------------------------------------------
# Multiple Signal Summary Statistics
#----------------------------------------------------------

# Only identify conditional summary statistics for loci with more than one signal

# Cis-pQTL
doParallel::registerDoParallel(cores=16)
foreach(locus=cs.cis.n.signals$Gene[cs.cis.n.signals$Signals > 1]) %dopar% {

    pqtl.locus <- cis.pqtl.loci[[locus]] %>%
        dplyr::mutate(N = NCHROBS / 2) %>%
        dplyr::select(SNP=SNP, A1=minor_allele, A2=major_allele, freq=MAF, b=pQTL_beta, pQTL_SE, p=pQTL_pval, N)

    write.table(pqtl.locus, paste0(locus, ".ma"), sep=" ", quote=F, row.names=F)

    # Save SNPs at pQTL locus
    write.table(pqtl.locus$SNP, paste0(locus, ".snps"), row.names=F, col.names=F, quote=F)

    # Save conditional SNP list for each conditional signal
    conditional.snps <- cs.cis.loci[[locus]]$SNP
    for (i in 1:length(conditional.snps)) {
        snp = conditional.snps[i]
        number = cs.cis.loci[[locus]]$Credibility_Set[i]
        write.table(setdiff(conditional.snps, snp), paste0(locus, "-", number, ".cond.snps"), row.names=F, col.names=F, quote=F)
    }
}

# Trans-pQTL
doParallel::registerDoParallel(cores=16)
foreach(locus=cs.trans.n.signals$Gene[cs.trans.n.signals$Signals > 1]) %dopar% {

    pqtl.locus <- trans.pqtl.loci[[locus]] %>%
        dplyr::mutate(N = NCHROBS / 2) %>%
        dplyr::select(SNP=SNP, A1=minor_allele, A2=major_allele, freq=MAF, b=pQTL_beta, pQTL_SE, p=pQTL_pval, N)

    write.table(pqtl.locus, paste0("Trans-", locus, ".ma"), sep=" ", quote=F, row.names=F)

    # Save SNPs at pQTL locus
    write.table(pqtl.locus$SNP, paste0("Trans-", locus, ".snps"), row.names=F, col.names=F, quote=F)

    # Save conditional SNP list for each conditional signal
    conditional.snps <- cs.trans.loci[[locus]]$SNP
    for (i in 1:length(conditional.snps)) {
        snp = conditional.snps[i]
        number = cs.trans.loci[[locus]]$Credibility_Set[i]
        write.table(setdiff(conditional.snps, snp), paste0("Trans-", locus, "-", number, ".cond.snps"), row.names=F, col.names=F, quote=F)
    }
}
