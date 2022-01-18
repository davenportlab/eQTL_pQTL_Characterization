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
library(coloc)
library(susieR)
library(parallel)

options(stringsAsFactors = FALSE)

args = commandArgs(trailingOnly = TRUE)
locus.name = args[1]
locus.input = args[2]
output = args[3]

#----------------------------------------------------------
# Load Data
#----------------------------------------------------------

eqtl.locus <- readRDS(locus.input)

#----------------------------------------------------------
# Identify Credible Sets
#----------------------------------------------------------

eqtl.S <- runsusie(
    eqtl.locus,
    prior_variance=(0.15/eqtl.locus$sdY)^2,
    estimate_prior_variance=T
)

credible.sets <- summary(eqtl.S)$vars %>%
    dplyr::filter(cs > 0) %>%
    dplyr::arrange(cs) %>%
    dplyr::mutate(snp = eqtl.locus$snp[variable]) %>%
    dplyr::mutate(gene = locus.name) %>%
    dplyr::select(gene, snp, snp_prob=variable_prob, credible_set=cs)

write.table(credible.sets, output, quote=F, sep="\t", row.names=F, col.names=F)
