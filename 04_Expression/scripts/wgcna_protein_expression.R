#----------------------------------------------------------
# WGCNA Analysis of Proteins
# Created: 14 October 2021
#----------------------------------------------------------

#----------------------------------------------------------
# Environment
#----------------------------------------------------------

rm(list=ls())

library(WGCNA)
library(tidyverse)
library(RColorBrewer)
library(biomaRt)

options(stringsAsFactors = FALSE)

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

protein.exp <- read.csv("04_Expression/data/protein_expression/data_291x1860_MS2019.csv", row.names=1)
protein.info <- read.csv("04_Expression/data/protein_expression/protein_info_291_MS2019.csv")
sample.info <- read.csv("04_Expression/data/protein_expression/sample_info_1860_MS2019.csv", row.names=1)

colnames(protein.exp) <- sapply(strsplit(colnames(protein.exp), "\\."), function(x) substring(x[1], 2))
protein.exp <- as.data.frame(t(protein.exp))

#----------------------------------------------------------
# WGCNA Checks
#----------------------------------------------------------

# Identifies excessive missing values and outliers
# Generally more important for microarray data
gsg <- goodSamplesGenes(protein.exp)
gsg$allOK

#----------------------------------------------------------
# Soft Threshold
#----------------------------------------------------------

# Identify the soft threshold power that will be used to generate the network

if (!file.exists("04_Expression/data/protein_expression_generated/TOM.matrix.RDS")) {
  
  powers = seq(1, 20, by=1)
  
  soft.thresholds = pickSoftThreshold(protein.exp, powerVector=powers, networkType="signed")
  
  soft.threshold = soft.thresholds$powerEstimate
}

#----------------------------------------------------------
# Adjacency and TOM Network
#----------------------------------------------------------

if (!file.exists("04_Expression/data/protein_expression_generated/TOM.matrix.RDS")) {
  
  adjacency.matrix = adjacency(protein.exp, type="signed", power=soft.threshold)
  
  TOM.matrix = TOMsimilarity(adjacency.matrix, TOMType="signed")
  
  saveRDS(TOM.matrix, file="04_Expression/data/protein_expression_generated/TOM.matrix.RDS")
  
} else {
  
  TOM.matrix <- readRDS("04_Expression/data/protein_expression_generated/TOM.matrix.RDS")
}

TOM.dist = 1 - TOM.matrix
rm(TOM.matrix)

#----------------------------------------------------------
# Clustering TOM Matrix
#----------------------------------------------------------

dendrogram = hclust(as.dist(TOM.dist), method="average")

dynamic.mods = lapply(0:4, function(x) {
  cutreeDynamic(
    dendro=dendrogram, distM=TOM.dist, deepSplit=x, pamRespectsDendro=FALSE,
    minClusterSize = 10
  )
})
rm(TOM.dist)

names(dynamic.mods) <- paste("Deep Split", 0:4)

# Generates clusters of the following sizes
sapply(dynamic.mods, function(x) length(table(x)))

#----------------------------------------------------------
# Distribution of Co-Expression Module Sizes
#----------------------------------------------------------

plot.data <- lapply(1:length(dynamic.mods), function(i, x, n) {
  as.data.frame(table(x[i])) %>%
    dplyr::mutate(Clusters=n[i])
}, x=dynamic.mods, n=sapply(dynamic.mods, function(x) length(table(x)))) %>%
  do.call(rbind, .) %>%
  dplyr::select(Module=1, Frequency=Freq, Clusters)

ggplot(plot.data) +
  geom_col(aes(x=Module, y=Frequency), width=1) +
  xlab("Module") + ylab("Number of Genes") +
  facet_grid(Clusters~.) +
  ggplot_theme +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())
