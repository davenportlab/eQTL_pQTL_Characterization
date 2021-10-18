#----------------------------------------------------------
# Co-Expression Module Annotation
# Created: 16 October 2021
#----------------------------------------------------------

#----------------------------------------------------------
# Environment
#----------------------------------------------------------

rm(list=ls())

library(tidyverse)
library(enrichR)

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

gene.info <- read.table("04_Expression/data/gene_expression/Gene_info_864_20416.txt")
modules <- read.csv("04_Expression/data/gene_expression_generated/modules.csv")
eigengenes <- read.csv("04_Expression/data/gene_expression_generated/eigengenes.csv", row.names=1)
variance.explained <- read.csv("04_Expression/data/gene_expression_generated/variance.explained.csv")

#----------------------------------------------------------
# Ontology/Pathway Enrichment
#----------------------------------------------------------

setEnrichrSite("Enrichr")

enrichr.dbs <- listEnrichrDbs()

enrichr.dbs[grepl("GO", enrichr.dbs$libraryName),]
enrichr.dbs[grepl("KEGG", enrichr.dbs$libraryName),]
enrichr.dbs[grepl("Reactome", enrichr.dbs$libraryName),]

selected.dbs <- c(
  "GO_Biological_Process_2021",
  "GO_Cellular_Component_2021",
  "GO_Molecular_Function_2021",
  "KEGG_2021_Human",
  "Reactome_2016"
)

module.names <- paste0("Module_", 1:dim(eigengenes)[2])
module.list <- lapply(module.names, function(module.name) {
  modules %>%
    dplyr::filter(Module==module.name) %>%
    merge(., gene.info, by.x="Gene", by.y="gene_id") %>%
    dplyr::select(Gene.ID=Gene, Gene.Name=gene_name)
})
names(module.list) <- module.names

module.annotations <- lapply(module.list, function(module) {

  enriched <- enrichr(module$Gene.Name, selected.dbs)
  
  annots <- lapply(selected.dbs, function(db) {
    as.data.frame(enriched[[db]])
  })
  names(annots) <- selected.dbs
  
  return(annots)
})

annotations.by.db <- lapply(selected.dbs, function(db) {
  
  annots <- lapply(1:length(module.names), function(i, m, n) {
    m[[i]][[db]] %>%
      dplyr::mutate(Module=n[i]) %>%
      dplyr::select(Module, everything())
  }, m=module.annotations, n=module.names)
  names(annots) <- module.names
  
  return(annots)
})
names(annotations.by.db) <- selected.dbs

#----------------------------------------------------------
# Save Annotations
#----------------------------------------------------------

for (db in selected.dbs) {
  
  output.annot <- do.call(rbind, annotations.by.db[[db]]) %>%
    dplyr::filter(Adjusted.P.value < 0.05)
  write.csv(output.annot, paste0("04_Expression/data/gene_coexpression_module_annotations/", db, ".csv"), row.names=F)
}
