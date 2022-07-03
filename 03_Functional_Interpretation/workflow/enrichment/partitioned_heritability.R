rm(list=ls())

library(tidyverse)
library(lme4qtl)
library(parallel)
library(Matrix)

# Format:
#   partitioned_heritability.R <Annotation>
args = commandArgs(trailingOnly=TRUE)
annotation = args[1]

# Reads GRMs made using GCTA, which are written in a binary format
# Code is provided within GCTA documentation
# Modified to create a symmetric semi-definite matrix from the lower triangular matrix
ReadGRMBin=function(prefix, AllN=F, size=4){
    sum_i=function(i){
        return(sum(1:i))
    }
    BinFileName=paste(prefix,".grm.bin",sep="")
    NFileName=paste(prefix,".grm.N.bin",sep="")
    IDFileName=paste(prefix,".grm.id",sep="")
    id = read.table(IDFileName)
    n=dim(id)[1]
    BinFile=file(BinFileName, "rb");
    grm=readBin(BinFile, n=n*(n+1)/2, what=numeric(0), size=size)
    NFile=file(NFileName, "rb");
        if(AllN==T){
        N=readBin(NFile, n=n*(n+1)/2, what=numeric(0), size=size)
        }
    else N=readBin(NFile, n=1, what=numeric(0), size=size)
    i=sapply(1:n, sum_i)
    close(BinFile)
    close(NFile)
  
    GRM <- matrix(0.0, nrow=n, ncol=n, dimnames=list(id$V1, id$V2))
    GRM[!lower.tri(GRM,diag=T)] <- grm[-i]
    GRM <- GRM + t(GRM)
    diag(GRM) <- grm[i]
    return(GRM)
}

# Load data used for eigengene mapping
map.data <- read.csv("/nfs/users/nfs_n/nm18/gains_team282/nikhil/expression/eigengene_sva/mapping_data.csv") %>%
    dplyr::mutate(GAinS.ID.Annotation=GAinS.ID) %>%
    dplyr::mutate(GAinS.ID.Other=GAinS.ID)

# Load GRM for annotation and non-annotation SNPs
annotation.grm <- ReadGRMBin(prefix=paste0("/nfs/users/nfs_n/nm18/gains_team282/epigenetics/enrichment/heritability/", annotation, "/annotation_snps_grm"))
annotation.grm <- pmin(pmax(annotation.grm, 0), 1)

other.grm <- ReadGRMBin(prefix=paste0("/nfs/users/nfs_n/nm18/gains_team282/epigenetics/enrichment/heritability/", annotation, "/other_snps_grm"))
other.grm <- pmin(pmax(other.grm, 0), 1)

# Covariates include: 7 Genotyping PCs, 20 PEER factors, Diagnosis, SRS1, Neutrophils, Monocytes, Lymphocytes
f <- paste0(
    c(paste0("PC", 1:7),
    paste0("PEER_", 1:20),
    "Diagnosis", "SRS1",
    "Neutrophils", "Monocytes", "Lymphocytes"
), collapse="+")

# Use first module eigengene from each module
me.set <- colnames(map.data)[grepl("^ME_[0-9]*_1", colnames(map.data))]

var.comp <- mclapply(me.set, function(me) {

    tryCatch({

        # Use lme4qtl to model relationship matrices as separate variance components
        m = relmatLmer(
            as.formula(paste0(me, " ~ ", f, " + (1|GAinS.ID) + (1|GAinS.ID.Annotation) + (1|GAinS.ID.Other)")),
            map.data,
            relmat=list(GAinS.ID.Annotation=annotation.grm, GAinS.ID.Other=other.grm)
        )

        # Extract variance explained by each component
        var.props = lme4qtl::VarProp(m) %>%
            as.data.frame() %>%
            dplyr::mutate(Eigengene=me, Annotation=annotation) %>%
            dplyr::select(Eigengene, Annotation, Component=grp, Variance=vcov, Proportion=prop) %>%
            dplyr::mutate(Component = recode(Component, GAinS.ID="Individual", GAinS.ID.Annotation="Annotation", GAinS.ID.Other="Other"))
        
        return(var.props)

    }, error = function(e) {

        var.props = data.frame(
            Eigengene=me, Annotation=annotation, Component=NA, Variance=NA, Proportion=NA
        )

        return(var.props)
    })
}, mc.cores=16) %>%
    do.call(rbind, .)

write.table(var.comp, paste0(annotation, ".variance_components.csv"), quote=F, sep=",", row.names=F, col.names=F)
