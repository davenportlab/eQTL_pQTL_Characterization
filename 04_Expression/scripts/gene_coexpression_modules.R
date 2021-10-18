#----------------------------------------------------------
# Gene Co-Expression Modules
# Created: 14 October 2021
#----------------------------------------------------------

#----------------------------------------------------------
# Environment
#----------------------------------------------------------

rm(list=ls())

library(tidyverse)
library(RColorBrewer)
library(cowplot)

ggplot_theme <- theme_bw(base_size=18) +
  theme(
    panel.border=element_blank(),
    panel.grid=element_blank(),
    axis.line.x.bottom=element_line(color="black", size=0.25),
    axis.line.y.left=element_line(color="black", size=0.25),
    legend.position="bottom",
    strip.background=element_rect(fill="#EEEEEE", color="white", size=0.25)
  )

#----------------------------------------------------------
# Load Data
#----------------------------------------------------------

modules <- read.csv("04_Expression/data/gene_expression_generated/modules.csv")
eigengenes <- read.csv("04_Expression/data/gene_expression_generated/eigengenes.csv", row.names=1)
variance.explained <- read.csv("04_Expression/data/gene_expression_generated/variance.explained.csv")

sample.info <- read.table("04_Expression/data/gene_expression/Sample_info_864.txt")
sample.info <- sample.info %>% 
  dplyr::filter(supplier_name %in% rownames(eigengenes))
rownames(sample.info) <- sample.info$supplier_name
sample.info <- sample.info[rownames(eigengenes),]

srs.info <- read.table("04_Expression/data/gene_expression/full-gains-SRS-predictions_mNN-RF.tsv", header=T)
srs.info <- srs.info %>%
  dplyr::filter(Assay=="RNA-seq")
row.names(srs.info) <- srs.info$Sample_id
srs.info <- srs.info[rownames(eigengenes),]

#----------------------------------------------------------
# Association with SRSq
#----------------------------------------------------------

# Exploratory Correlation Matrix

variables <- cbind(eigengenes, srs.info$SRSq)
colnames(variables)[ncol(variables)] <- "SRSq"

rho.values <- abs(cor(variables, method="spearman"))

rho.values.hclust <- hclust(as.dist(1 - rho.values))

rho.values %>%
  as.data.frame() %>%
  dplyr::mutate(Variable.2=rownames(.)) %>%
  tidyr::gather("Variable.1", "Rho", -Variable.2) %>%
  dplyr::mutate(Variable.1=factor(Variable.1, levels=colnames(rho.values)[rho.values.hclust$order])) %>%
  dplyr::mutate(Variable.2=factor(Variable.2, levels=colnames(rho.values)[rho.values.hclust$order])) %>%
  ggplot() +
  geom_raster(aes(x=Variable.1, y=Variable.2, fill=Rho)) +
  scale_fill_distiller(palette="Blues", limits=c(0, 1), direction=1) +
  ggplot_theme +
  theme(axis.title=element_blank(), axis.text.x=element_text(angle=270), legend.position="right")
ggsave("04_Expression/results/eigengene_association_with_srsq.svg", width=15, height=15)

# Perform Correlation Test with Spearman's Rho
#   Filtered by Benjamini-Hochberg P-Value < 0.05 and Rho >= 0.8

estimates <- lapply(colnames(eigengenes), function(eigengene.name) {
  test.result = cor.test(srs.info$SRSq, eigengenes[,eigengene.name], alternative="two.sided", method="spearman", exact=FALSE)
  return(c(test.result$estimate, test.result$p.value))
}) %>%
  do.call(rbind, .) %>%
  as.data.frame() %>%
  dplyr::select(Rho=1, P.Value=2) %>%
  dplyr::mutate(Eigengene=colnames(eigengenes)) %>%
  dplyr::mutate(Adjusted.P.Value=p.adjust(P.Value, method="BH")) %>%
  dplyr::filter(abs(Rho) >= 0.8, Adjusted.P.Value < 0.05) %>%
  dplyr::arrange(Rho) %>%
  dplyr::mutate(Eigengene=factor(Eigengene, levels=Eigengene))
print(estimates)

# Plot Associated Eigengenes by SRS Group

plots <- lapply(estimates$Eigengene, function(eigengene.name) {
  cbind(srs.info$SRS, eigengenes) %>%
    as.data.frame() %>%
    dplyr::select(SRS=1, everything()) %>%
    ggplot(aes_string(x="SRS", y=as.character(eigengene.name))) +
    geom_boxplot() +
    ggplot_theme +
    theme(axis.title.x=element_blank(), axis.text.x=element_text(angle=30, hjust=1))
})

plot_grid(plotlist=plots)
ggsave("04_Expression/results/eigengene_association_with_srs_group.svg", width=12, height=10)
