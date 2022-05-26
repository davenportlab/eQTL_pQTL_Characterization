library(tidyverse)
library(parallel)

cis.overlaps <- read.table("cis_overlaps.tsv", header=T, row.names=1, check.names=F)
gains.overlaps <- read.table("gains_overlaps.tsv", header=T, row.names=1, check.names=F)
shared.overlaps <- read.table("shared_overlaps.tsv", header=T, row.names=1, check.names=F)

#----------------------------------------------------------
# Fisher Exact Tests
#----------------------------------------------------------

n.gains.snps <- nrow(gains.overlaps)
n.shared.snps <- nrow(shared.overlaps)

n.gains.overlaps <- colSums(gains.overlaps)
n.shared.overlaps <- colSums(shared.overlaps)

fisher.results <- do.call(rbind, lapply(1:ncol(gains.overlaps), function(i) {
    test = fisher.test(
        matrix(c(
            n.gains.overlaps[i], n.gains.snps - n.gains.overlaps[i], 
            n.shared.overlaps[i], n.shared.snps - n.shared.overlaps[i]
        ), nrow=2)
    )
    data.frame(P_Value=test$p.value, Odds_Ratio=test$estimate, Odds_Ratio_Lower=test$conf.int[1], Odds_Ratio_Upper=test$conf.int[2])
})) %>%
    dplyr::arrange(desc(Odds_Ratio))
rownames(fisher.results) <- colnames(gains.overlaps)

write.csv(fisher.results, "fisher_results.csv", quote=F)

#----------------------------------------------------------
# Controlling for Distance to TSS
#----------------------------------------------------------

cis.snps <- read.table("cis_esnps.sorted.bed")
colnames(cis.snps) <- c("chr", "start", "end", "ID", "dist_to_tss")
gains.snps <- read.table("gains_esnps.sorted.bed")
colnames(gains.snps) <- c("chr", "start", "end", "ID", "dist_to_tss")
shared.snps <- read.table("shared_esnps.sorted.bed")
colnames(shared.snps) <- c("chr", "start", "end", "ID", "dist_to_tss")

cis.snps <- cis.snps %>% 
    dplyr::group_by(ID) %>%
    dplyr::summarize(mean_dist_to_tss=mean(dist_to_tss)) %>%
    dplyr::mutate(quantile=cut(mean_dist_to_tss, breaks=seq(0, 10^6, length.out=50), include_lowest=TRUE, right=FALSE))

gains.snps <- gains.snps %>% 
    dplyr::group_by(ID) %>%
    dplyr::summarize(mean_dist_to_tss=mean(dist_to_tss)) %>%
    dplyr::mutate(quantile=cut(mean_dist_to_tss, breaks=seq(0, 10^6, length.out=50), include_lowest=TRUE, right=FALSE))

shared.snps <- shared.snps %>% 
    dplyr::group_by(ID) %>%
    dplyr::summarize(mean_dist_to_tss=mean(dist_to_tss)) %>%
    dplyr::mutate(quantile=cut(mean_dist_to_tss, breaks=seq(0, 10^6, length.out=50), include_lowest=TRUE, right=FALSE))

gains.snp.tss.freq <- table(gains.snps$quantile)
shared.snp.tss.freq <- table(shared.snps$quantile)

set.seed(42)
overlap.samples <- do.call(rbind, mclapply(1:10000, function(i) {
    snp.sample <- do.call(rbind, lapply(
        names(gains.snp.tss.freq),
        function(x) {
            cis.snps %>%
                dplyr::filter(quantile == x) %>%
                dplyr::sample_n(size=as.numeric(gains.snp.tss.freq[x]), replace=TRUE)
        }
    ))
    colSums(cis.overlaps[snp.sample$ID, ])
}, mc.cores=16))

gains.tss.corrected.results <- data.frame(
    rowSums(apply(overlap.samples, 1, function(x) { x >= n.gains.overlaps })) / 10000
) %>%
    dplyr::select(P_Value=1) %>%
    dplyr::arrange(P_Value)

write.csv(gains.tss.corrected.results, "gains_tss_corrected_results.csv", quote=F)

set.seed(42)
overlap.samples <- do.call(rbind, mclapply(1:10000, function(i) {
    snp.sample <- do.call(rbind, lapply(
        names(shared.snp.tss.freq),
        function(x) {
            cis.snps %>%
                dplyr::filter(quantile == x) %>%
                dplyr::sample_n(size=as.numeric(shared.snp.tss.freq[x]), replace=TRUE)
        }
    ))
    colSums(cis.overlaps[snp.sample$ID, ])
}, mc.cores=16))

shared.tss.corrected.results <- data.frame(
    rowSums(apply(overlap.samples, 1, function(x) { x >= n.shared.overlaps })) / 10000
) %>%
    dplyr::select(P_Value=1) %>%
    dplyr::arrange(P_Value)

write.csv(shared.tss.corrected.results, "shared_tss_corrected_results.csv", quote=F)
