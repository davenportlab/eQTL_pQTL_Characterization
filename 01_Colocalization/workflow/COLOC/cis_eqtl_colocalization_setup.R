library(tidyverse)
library(GenomicRanges)
library(data.table)

cis.eqtl <- read.table("/nfs/users/nfs_n/nm18/gains_team282/eqtl/cisresults/eigenMT/ciseqtl_eigenMT_corrected.txt") %>%
    dplyr::filter(Sig)

cis.eqtl.sig <- readRDS("/nfs/users/nfs_n/nm18/gains_team282/eqtl/cisresults/cisqtl_all_significant.rds")

gene.info <- read.table("/nfs/team282/data/gains_team282/gene_info_864_20412_hla.txt", sep="\t") %>%
    dplyr::filter(gene_id %in% cis.eqtl.sig$gene) %>%
    dplyr::mutate(TSS = ifelse(strand == "+", start, end)) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(TSS.Start = max(TSS - 1e6, 1)) %>%
    dplyr::mutate(TSS.End = TSS + 1e6 - 1) %>%
    dplyr::filter(!is.na(seqnames), !is.na(TSS), !is.na(TSS.Start), !is.na(TSS.End))

tss.ranges <- makeGRangesFromDataFrame(
    gene.info, 
    seqnames.field="seqnames",
    start.field="TSS.Start",
    end.field="TSS.End"
)

overlaps <- findOverlaps(tss.ranges, tss.ranges)
