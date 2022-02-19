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
library(lme4)

options(stringsAsFactors = FALSE)

args = commandArgs(trailingOnly = TRUE)
locus.name = args[1]
locus.input = args[2]
genotypes.input = args[3]
expression.input = args[4]
covariates.input = args[5]
output = args[6]

#----------------------------------------------------------
# Load Data
#----------------------------------------------------------

eqtl.locus <- readRDS(locus.input)
genotypes <- readRDS(genotypes.input)
expression <- readRDS(expression.input)
covariates <- readRDS(covariates.input)

#----------------------------------------------------------
# Identify Credible Sets
#----------------------------------------------------------

# Create a full design matrix with all genotypes
full.mtx <- cbind(expression, covariates, genotypes)

covs <- setdiff(colnames(covariates), c("GAinS.ID"))

summary.stats <- do.call(rbind, mclapply(colnames(genotypes), function(snp) {

    # Regress expression against fixed effects (SNP and covariates) and random effects (patient)
    # This is the same model used for eQTL mapping
    variant.design <- full.mtx[,c(colnames(expression), snp, colnames(covariates))]

    f.alt <- as.formula(paste0(locus.name, "~`", snp , "`+", paste0(covs, collapse="+"), "+(1|GAinS.ID)"))
    model.test <- lmer(f.alt, data=variant.design, REML=FALSE)

    if (!all(complete.cases(variant.design[, snp]))) {
        model.test <- update(model.test, subset=complete.cases(variant.design[, snp]))
    }

    # Regress expression against SNP, covariates, and intercept of random effects from the initial model
    regress.design <- merge(
        variant.design %>% dplyr::mutate(Sample.ID=rownames(.)), 
        coef(model.test)$GAinS.ID[,1,drop=F], 
        by.x="GAinS.ID", by.y=0
    )
    colnames(regress.design)[ncol(regress.design)] <- "GAinS.ID.Intercept"

    regress.design[,snp] = regress.design[,snp] - regress.design$GAinS.ID.Intercept

    f.regress <- as.formula(paste0(locus.name, "~`", snp, "`+", paste0(covs, collapse="+")))
    model.regress <- lm(f.regress, data=regress.design)

    return(summary(model.regress)$coefficients[snp, c("Estimate", "Std. Error")])
}))

LD <- cor(genotypes, use="pairwise.complete.obs")

z.scores <- summary.stats[,1] / summary.stats[,2]

tryCatch({

    susie.output <- susie_rss(z.scores, LD)
    
    credible.sets <- summary(susie.output)$vars %>%
        dplyr::filter(cs > 0) %>%
        dplyr::arrange(cs, desc(variable_prob)) %>%
        dplyr::mutate(snp = eqtl.locus$snp[variable]) %>%
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
