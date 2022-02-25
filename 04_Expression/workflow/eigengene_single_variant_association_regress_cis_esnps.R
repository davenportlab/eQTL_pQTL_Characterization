#----------------------------------------------------------
# Author: Nikhil Milind
# Created: 25 November 2021
#----------------------------------------------------------

#----------------------------------------------------------
# Environment
#----------------------------------------------------------

library(tidyverse)
library(data.table)

args = commandArgs(trailingOnly=TRUE)
genotypes.file <- args[1]
design.matrix.file <- args[2]

#----------------------------------------------------------
# Load mQTL Data
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
# Load Expression Data
#----------------------------------------------------------

# Load gene expression data
gene.exp <- read.table("/lustre/scratch119/humgen/projects/gains_team282/eqtl/data/logcpm_864_20412_hla.txt") %>%
    t() %>% as.data.frame()
rownames(gene.exp) <- gsub("^GA", "", rownames(gene.exp))
gene.exp <- gene.exp[design.matrix$Sample.ID, ]

# Regression design
regress.design <- cbind(gene.exp, genotypes)

# Load lead SNPs that need to be regressed out
lead.cis.eqtl <- read.table("/nfs/users/nfs_n/nm18/gains_team282/eqtl/cisresults/eigenMT/ciseqtl_eigenMT_corrected.txt") %>%
    dplyr::filter(Sig) %>%
    dplyr::filter(snps %in% colnames(genotypes))

# Regress out each cis-eSNP
for (i in 1:nrow(lead.cis.eqtl)) {
    regress.gene <- lead.cis.eqtl$gene[i]
    regress.snp <- lead.cis.eqtl$snps[i]
    lm.formula <- as.formula(paste0(regress.gene, "~`", regress.snp, "`"))
    lm.fit <- lm(lm.formula, data=regress.design)
    lm.resids <- regress.design[, regress.gene] - predict(lm.fit, newdata=regress.design, na.action=na.pass)
    regress.design[, regress.gene] <- lm.resids
}

# Load module membership information
modules <- read.csv("~/gains_team282/nikhil/expression/gene_expression/modules.csv")

genes.by.module <- split(modules$Gene, modules$Module)

eigengenes <- do.call(cbind, lapply(genes.by.module, function(genes) {
    
    module.exp = regress.design[, genes]
    gene.means = apply(module.exp, 2, function(x) { mean(x, na.rm=TRUE) })
    module.exp.imputed = t(apply(module.exp, 1, function(x) { x[is.na(x)] = gene.means[is.na(x)]; return(x) }))
    svd(scale(module.exp.imputed))$u[,1]
}))
colnames(eigengenes) <- gsub("Module_", "ME_", names(genes.by.module))
rownames(eigengenes) <- rownames(gene.exp)

# Save regressed eigengenes
saveRDS(eigengenes, "regressed.eigengenes.RDS")
