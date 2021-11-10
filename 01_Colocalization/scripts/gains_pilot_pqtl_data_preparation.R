#----------------------------------------------------------
# Prepare GAinS Pilot pQTL for Colocalization
# Created: 25 October 2021
#----------------------------------------------------------

#----------------------------------------------------------
# Environment
#----------------------------------------------------------

rm(list=ls())

library(tidyverse)

options(stringsAsFactors = FALSE)

setwd("~/eQTL_pQTL_Characterization/")

source("01_Colocalization/scripts/utils/ggplot_theme.R")

#----------------------------------------------------------
# Load pQTL Data
#----------------------------------------------------------

cis.pqtl <- readRDS("~/gains_team282/proteomics/pqtl/fp_pilot/cis_pqtl_fp49.rds")
cis.pqtl.proteins <- sapply(sapply(strsplit(as.character(cis.pqtl$Gene), "\\|"), function(x) strsplit(x[2], "-")), function(x) x[1])

uniprot.id.map <- read.table("01_Colocalization/data/UniProt/HUMAN_9606_idmapping_selected.tab", sep="\t", stringsAsFactors=F)

uniprot.id.map <- uniprot.id.map %>%
  dplyr::select(UniProt.ID=1, Ensembl.ID=19) %>%
  dplyr::filter(UniProt.ID %in% cis.pqtl.proteins) %>%
  dplyr::filter(Ensembl.ID != "")

# Some IDs have been retired or changed
old.accessions <- data.frame(
  UniProt.ID=c(
    "P18464", "A0A087WW79", "A0A075B6N7", "A0A087X188",
    "H0YLH7", "A0A087WXL1", "A0A087WU43", "A0A0A0MR60",
    "J3KNP4", "A0A087WW89", "A0A087WXW9", "A0A087X1N7",
    "A0A087WW55", "A0A075B6H6", "A0A0A0MQY0", "A0A0A0MTS2",
    "A0A0C4DFU2", "A0A0U1RR20", "P01860", "P0CG04",
    "A0M8Q6", "A8K878"
  ),
  Ensembl.ID=c(
    "ENSG00000234745", "ENSG00000006432", "ENSG00000211890", "ENSG00000110934",
    "ENSG00000140564", "ENSG00000110203", "ENSG00000039068", "ENSG00000105426",
    "ENSG00000185033", "ENSG00000225698", "ENSG00000130635", "ENSG00000183091",
    "ENSG00000204983", "ENSG00000211592", "ENSG00000154553", "ENSG00000105220",
    "ENSG00000112096", "ENSG00000116690", "ENSG00000211897", "ENSG00000211675",
    "ENSG00000211685", "ENSG00000145050"
  ),
  stringsAsFactors = FALSE
)

gene.exp <- read.table("04_Expression/data/gene_expression/Logcpm_864_20417_HLA.txt")

uniprot.id.map <- uniprot.id.map %>%
  dplyr::bind_rows(old.accessions) %>%
  dplyr::mutate(Ensembl.ID.in.Gene.Exp=sapply(lapply(strsplit(Ensembl.ID, "; "), function(x) x[x %in% rownames(gene.exp)]), function(x) paste0(x, collapse="; ")))

rownames(uniprot.id.map) <- uniprot.id.map$UniProt.ID

cis.pqtl$Ensembl.ID = uniprot.id.map[cis.pqtl.proteins, "Ensembl.ID.in.Gene.Exp"]

cis.pqtl <- cis.pqtl %>%
  dplyr::mutate(varbeta=eQTL_SE^2) %>%
  dplyr::mutate(SNP=as.character(SNP)) %>%
  dplyr::select(Gene, Ensembl.ID, snp=SNP, beta=eQTL_beta, varbeta)

cis.pqtl.loci <- split(cis.pqtl[,c("snp", "beta", "varbeta", "Ensembl.ID")], cis.pqtl$Gene)

saveRDS(cis.pqtl.loci, "~/gains_team282/nikhil/colocalization/pilot.cis.pqtl.loci.RDS")

#----------------------------------------------------------
# Load Gene Expression Data
#----------------------------------------------------------

protein.exp <- read.table("~/gains_team282/proteomics/pqtl/fp_pilot/log.norm.flt.imp.data.txt")

protein.exp.sd <- data.frame(
  Gene=rownames(protein.exp),
  SD=apply(protein.exp, 1, sd)
)

protein.exp.sd <- protein.exp.sd[names(cis.pqtl.loci),]

saveRDS(protein.exp.sd, "~/gains_team282/nikhil/colocalization/pilot.cis.pqtl.sdY.RDS")
