#----------------------------------------------------------
# Gene Co-Expression Modules
# Created: 14 October 2021
#----------------------------------------------------------

#----------------------------------------------------------
# Environment
#----------------------------------------------------------

rm(list=ls())

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
# Association with SRS Group
#----------------------------------------------------------

cors <- abs(cor(cbind(eigengenes, srs.info$SRSq)))
cors[,"srs.info$SRSq"][order(cors[,"srs.info$SRSq"], decreasing=TRUE)]

plot.data <- cbind(eigengenes, srs.info)
plot.data$SRS <- factor(plot.data$SRS)

ggplot(plot.data) + 
  geom_boxplot(aes(x=SRS, y=ME_7))
