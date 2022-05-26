library(tidyverse)
library(data.table)

cis.eqtl <- read.table("/nfs/users/nfs_n/nm18/gains_team282/eqtl/cisresults/eigenMT/ciseqtl_eigenMT_corrected.txt") %>%
    dplyr::filter(Sig)

categories  <- read.table("/lustre/scratch119/humgen/projects/gains_team282/eqtl/gtex/sepsis_specificity_categories.txt", row.names=1, sep="\t")

geno.bim <- fread("/lustre/scratch119/humgen/projects/gains_team282/Genotyping/All_genotyping_merged_filtered_b38_refiltered_rsID.bim") %>%
    as.data.frame()
colnames(geno.bim) <- c("chr", "snp", "cM", "pos", "minor", "major")

ld <- read.table("/nfs/users/nfs_n/nm18/gains_team282/nikhil/colocalization/cis_eqtl/fine_mapping/LD/lead_snps.80r2.tags.tsv", header=T)

tss.info <- read.table("/nfs/team282/data/gains_team282/Gene_info_864_20416.txt") %>%
    dplyr::mutate(tss = ifelse(strand == "+", start, end)) %>%
    dplyr::select(gene_id, tss)

# All cis-eQTL

cis.esnps <- cis.eqtl %>%
    merge(., ld, by.x="snps", by.y="SNP") %>%
    dplyr::mutate(TAGS=ifelse(TAGS == "NONE", snps, paste0(snps, "|", TAGS))) %>%
    dplyr::mutate(gene_snp_id = paste0(gene, "-", snps)) %>%
    dplyr::select(gene_snp_id, tags=TAGS)

cis.esnps.tags <- lapply(strsplit(cis.esnps$tags, "\\|"), function(x) { t(rbind(x)) })
names(cis.esnps.tags) <- cis.esnps$gene_snp_id

cis.esnps.bed <- do.call(rbind, 
    lapply(names(cis.esnps.tags), function(esnp.id) { cis.esnps.tags[[esnp.id]] %>% as.data.frame() %>% dplyr::mutate(ID=esnp.id) })
) %>%
    dplyr::select(SNP=1, ID) %>%
    dplyr::mutate(Gene=gsub("-.*", "", ID)) %>%
    merge(tss.info, ., by.x="gene_id", by.y="Gene") %>%
    merge(geno.bim, ., by.x="snp", by.y="SNP") %>%
    dplyr::mutate(start = pos - 1, end = pos, dist_to_tss=abs(tss - pos)) %>%
    dplyr::select(chr, start, end, ID, dist_to_tss)

fwrite(cis.esnps.bed, "cis_esnps.bed", sep="\t", col.names=F, scipen=99)

# Shared eSNPs

shared.esnps <- categories %>% 
    dplyr::filter(shared == "shared") %>%
    merge(., ld, by.x="snps", by.y="SNP") %>%
    dplyr::mutate(TAGS=ifelse(TAGS == "NONE", snps, paste0(snps, "|", TAGS))) %>%
    dplyr::mutate(gene_snp_id = paste0(gene, "-", snps)) %>%
    dplyr::select(gene_snp_id, tags=TAGS)

shared.esnps.tags <- lapply(strsplit(shared.esnps$tags, "\\|"), function(x) { t(rbind(x)) })
names(shared.esnps.tags) <- shared.esnps$gene_snp_id

shared.esnps.bed <- do.call(rbind, 
    lapply(names(shared.esnps.tags), function(esnp.id) { shared.esnps.tags[[esnp.id]] %>% as.data.frame() %>% dplyr::mutate(ID=esnp.id) })
) %>%
    dplyr::select(SNP=1, ID) %>%
    dplyr::mutate(Gene=gsub("-.*", "", ID)) %>%
    merge(tss.info, ., by.x="gene_id", by.y="Gene") %>%
    merge(geno.bim, ., by.x="snp", by.y="SNP") %>%
    dplyr::mutate(start = pos - 1, end = pos, dist_to_tss=abs(tss - pos)) %>%
    dplyr::select(chr, start, end, ID, dist_to_tss)

fwrite(shared.esnps.bed, "shared_esnps.bed", sep="\t", col.names=F, scipen=99)

# GAinS eSNPs

gains.esnps <- categories %>% 
    dplyr::filter(shared == "Bigger effect in GAinS") %>%
    merge(., ld, by.x="snps", by.y="SNP") %>%
    dplyr::mutate(TAGS=ifelse(TAGS == "NONE", snps, paste0(snps, "|", TAGS))) %>%
    dplyr::mutate(gene_snp_id = paste0(gene, "-", snps)) %>%
    dplyr::select(gene_snp_id, tags=TAGS)

gains.esnps.tags <- lapply(strsplit(gains.esnps$tags, "\\|"), function(x) { t(rbind(x)) })
names(gains.esnps.tags) <- gains.esnps$gene_snp_id

gains.esnps.bed <- do.call(rbind, 
    lapply(names(gains.esnps.tags), function(esnp.id) { gains.esnps.tags[[esnp.id]] %>% as.data.frame() %>% dplyr::mutate(ID=esnp.id) })
) %>%
    dplyr::select(SNP=1, ID) %>%
    dplyr::mutate(Gene=gsub("-.*", "", ID)) %>%
    merge(tss.info, ., by.x="gene_id", by.y="Gene") %>%
    merge(geno.bim, ., by.x="snp", by.y="SNP") %>%
    dplyr::mutate(start = pos - 1, end = pos, dist_to_tss=abs(tss - pos)) %>%
    dplyr::select(chr, start, end, ID, dist_to_tss)

fwrite(gains.esnps.bed, "gains_esnps.bed", sep="\t", col.names=F, scipen=99)
