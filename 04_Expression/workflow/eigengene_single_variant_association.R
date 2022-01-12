#----------------------------------------------------------
# Author: Nikhil Milind
# Created: 25 November 2021
#----------------------------------------------------------

#----------------------------------------------------------
# Environment
#----------------------------------------------------------

if (!require("lme4")) {
    install.packages("lme4")
}

library(data.table)
library(lme4)
library(parallel)

args = commandArgs(trailingOnly=TRUE)
genotypes.file <- args[1]
design.matrix.file <- args[2]
ME <- args[3]
output <- args[4]

#----------------------------------------------------------
# Load Data
#----------------------------------------------------------

# Design Matrix
design.matrix <- read.csv(design.matrix.file)

# Genotype Matrix
genotypes <- fread(genotypes.file, sep=" ", drop=2:6)

# Clean Genotype Matrix
patient.sample.match <- match(design.matrix$GAinS.ID, genotypes$FID)
genotypes <- genotypes[patient.sample.match,]
colnames(genotypes) <- gsub("X", "", colnames(genotypes))
colnames(genotypes) <- sapply(strsplit(colnames(genotypes), "_"), function(x) x[1])
genotypes[, 1] <- NULL
genotypes <- as.matrix(genotypes)
rownames(genotypes) <- design.matrix$Sample.ID

genotype.ids <- colnames(genotypes)

#----------------------------------------------------------
# Linear Mixed Model
#----------------------------------------------------------

# Create a full design matrix with all genotypes
genotypes <- cbind(design.matrix, genotypes)

all.vars <- colnames(design.matrix)
eigens <- all.vars[grepl("ME", all.vars)]
covs <- setdiff(setdiff(all.vars, eigens), c("Sample.ID", "GAinS.ID"))

results <- rbindlist(mclapply(genotype.ids, function(snp) {

    variant.design <- genotypes[,c(ME, snp, covs, "GAinS.ID")]

    f.null <- as.formula(paste0(ME, "~", paste0(covs, collapse="+"), "+(1|GAinS.ID)"))
    model.null <- lmer(f.null, data=variant.design, REML=FALSE)

    f.alt <- as.formula(paste0(ME, "~`", snp , "`+", paste0(covs, collapse="+"), "+(1|GAinS.ID)"))
    model.test <- lmer(f.alt, data=variant.design, REML=FALSE)

    if (!all(complete.cases(variant.design[, snp]))) {
        model.null <- update(model.null, subset=complete.cases(variant.design[, snp]))
        model.test <- update(model.test, subset=complete.cases(variant.design[, snp]))
    }

    data.frame(matrix(
        data=c(
            summary(model.test)$coefficients[2, ],  # Second row is for the SNP
            anova(model.null, model.test)["model.test", "Pr(>Chisq)"]
        ),
        nrow=1, ncol=4
    ))
}))

#----------------------------------------------------------
# Write Output
#----------------------------------------------------------

# Results has the following columns:
#   SNP
#   Beta
#   Standard Error
#   t Value
#   P-Value from ANOVA

results <- cbind(genotype.ids, results)

write.table(results, output, quote=F, row.names=F, col.names=F, sep="\t")
