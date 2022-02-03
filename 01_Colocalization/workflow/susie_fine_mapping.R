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
library(parallel)

options(stringsAsFactors = FALSE)

args = commandArgs(trailingOnly = TRUE)
locus.name = args[1]
locus.input = args[2]
genotypes.input = args[3]
var_y.input = args[4]
n_samples.input = args[5]
output = args[6]

#----------------------------------------------------------
# Load Data
#----------------------------------------------------------

eqtl.locus <- readRDS(locus.input)
genotypes <- readRDS(genotypes.input)
var.y <- readRDS(var_y.input)
n.samples <- readRDS(n_samples.input)

#----------------------------------------------------------
# Identify Credible Sets
#----------------------------------------------------------

LD <- cor(genotypes, use="pairwise.complete.obs")

tryCatch({

    susie.output <- susie_suff_stat(
        bhat=eqtl.locus$beta,
        shat=eqtl.locus$se,
        R=LD,
        n=n.samples,
        var_y=var.y
    )
    
    credible.sets <- summary(susie.output)$vars %>%
        dplyr::filter(cs > 0) %>%
        dplyr::arrange(cs) %>%
        dplyr::mutate(snp = eqtl.locus$snp[variable]) %>%
        dplyr::mutate(gene = locus.name) %>%
        dplyr::mutate(notes = "") %>%
        dplyr::select(gene, snp, snp_prob=variable_prob, credible_set=cs, notes)

    write.table(credible.sets, output, quote=F, sep="\t", row.names=F, col.names=F)

}, error = function(e) {

    credible.sets <- data.frame(
        gene=locus.name, snp=NA, snp_prob=NA, credible_set=NA, notes="SuSiE Error"
    )

    write.table(credible.sets, output, quote=F, sep="\t", row.names=F, col.names=F)
})
