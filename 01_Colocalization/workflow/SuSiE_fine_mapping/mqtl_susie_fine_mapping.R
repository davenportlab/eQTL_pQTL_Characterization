#----------------------------------------------------------
# Fine Mapping with SuSiE
# Created: 14 December 2021
#----------------------------------------------------------

#----------------------------------------------------------
# Environment
#----------------------------------------------------------

rm(list=ls())

library(tidyverse)
library(data.table)
library(susieR)

options(stringsAsFactors = FALSE)

args = commandArgs(trailingOnly = TRUE)
locus.name = args[1]
locus.input = args[2]
genotypes.input = args[3]
var.y.input = args[4]
n.samples.input = args[5]
signals = as.numeric(args[6])
output = args[7]

#----------------------------------------------------------
# Load Data
#----------------------------------------------------------

mqtl.locus <- readRDS(locus.input)
genotypes <- readRDS(genotypes.input)
var.y <- readRDS(var.y.input)
n.samples <- readRDS(n.samples.input)

#----------------------------------------------------------
# Identify Credible Sets
#----------------------------------------------------------

LD <- cor(genotypes, use="pairwise.complete.obs")

tryCatch({

    susie.output <- susie_suff_stat(
        bhat=mqtl.locus$Beta,
        shat=mqtl.locus$SE,
        R=LD,
        n=n.samples,
        var_y=var.y,
        L=signals
    )
    
    credible.sets <- summary(susie.output)$vars %>%
        dplyr::filter(cs > 0) %>%
        dplyr::arrange(cs, desc(variable_prob)) %>%
        dplyr::mutate(snp = mqtl.locus$SNP[variable]) %>%
        dplyr::mutate(gene = locus.name) %>%
        dplyr::mutate(notes = rep(NA, nrow(.))) %>%
        dplyr::select(gene, snp, snp_prob=variable_prob, credible_set=cs, notes)

    if (nrow(credible.sets) == 0) {
        credible.sets <- data.frame(
            gene=locus.name, snp=NA, snp_prob=NA, credible_set=NA, notes="No Credibility Sets"
        )
    }

    write.table(credible.sets, output, quote=F, sep="\t", row.names=F, col.names=F)

}, error = function(e) {

    credible.sets <- data.frame(
        gene=locus.name, snp=NA, snp_prob=NA, credible_set=NA, notes="SuSiE Error"
    )

    write.table(credible.sets, output, quote=F, sep="\t", row.names=F, col.names=F)
})
