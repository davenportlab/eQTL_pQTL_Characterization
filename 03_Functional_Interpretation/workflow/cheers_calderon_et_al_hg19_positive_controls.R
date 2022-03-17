library(tidyverse)

# Fairfax et al. 2012 - Nature Genetics (DOI: 10.1038/ng.2205)

mono.eqtl.info <- readxl::read_xls("~/eQTL_pQTL_Characterization/03_Functional_Interpretation/data/41588_2012_BFng2205_MOESM2_ESM.xls", sheet=2)
b.cell.eqtl.info <- readxl::read_xls("~/eQTL_pQTL_Characterization/03_Functional_Interpretation/data/41588_2012_BFng2205_MOESM2_ESM.xls", sheet=3)

mono.eSNPs <- mono.eqtl.info %>%
    dplyr::group_by(Gene) %>%
    dplyr::slice_min(order_by=P, n=1) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(CHR = paste0("chr", CHR)) %>%
    dplyr::select(SNP, CHR, BP) %>%
    unique()

b.cell.eSNPs <- b.cell.eqtl.info %>%
    dplyr::group_by(Gene) %>%
    dplyr::slice_min(order_by=P, n=1) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(CHR = paste0("chr", CHR)) %>%
    dplyr::select(SNP, CHR, BP) %>%
    unique()

write.table(mono.eSNPs, "snps_lists/monocyte.txt", quote=F, col.names=F, row.names=F, sep="\t")
write.table(b.cell.eSNPs, "snps_lists/b_cell.txt", quote=F, col.names=F, row.names=F, sep="\t")

# Fairfax et al. 2014 - Science (DOI: 10.1126/science.1246949)

mono.sim.eqtl.info <- readxl::read_xlsx("~/eQTL_pQTL_Characterization/03_Functional_Interpretation/data/1246949stables2.xlsx", col_types="text", sheet=2)

mono.lps.2.eSNPs <- mono.sim.eqtl.info %>%
    dplyr::select(SNP, Gene, SNP.Chrm, SNP.pos, LPS2.p.value) %>%
    dplyr::group_by(Gene) %>%
    dplyr::slice_min(order_by=LPS2.p.value) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(SNP.Chrm = paste0("chr", SNP.Chrm)) %>%
    dplyr::select(SNP, SNP.Chrm, SNP.pos) %>%
    unique()

mono.lps.24.eSNPs <- mono.sim.eqtl.info %>%
    dplyr::select(SNP, Gene, SNP.Chrm, SNP.pos, LPS24.p.value) %>%
    dplyr::group_by(Gene) %>%
    dplyr::slice_min(order_by=LPS24.p.value) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(SNP.Chrm = paste0("chr", SNP.Chrm)) %>%
    dplyr::select(SNP, SNP.Chrm, SNP.pos) %>%
    unique()

mono.ifn.eSNPs <- mono.sim.eqtl.info %>%
    dplyr::select(SNP, Gene, SNP.Chrm, SNP.pos, IFN.p.value) %>%
    dplyr::group_by(Gene) %>%
    dplyr::slice_min(order_by=IFN.p.value) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(SNP.Chrm = paste0("chr", SNP.Chrm)) %>%
    dplyr::select(SNP, SNP.Chrm, SNP.pos) %>%
    unique()

mono.naive.eSNPs <- mono.sim.eqtl.info %>%
    dplyr::select(SNP, Gene, SNP.Chrm, SNP.pos, Naive.p.value) %>%
    dplyr::group_by(Gene) %>%
    dplyr::slice_min(order_by=Naive.p.value) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(SNP.Chrm = paste0("chr", SNP.Chrm)) %>%
    dplyr::select(SNP, SNP.Chrm, SNP.pos) %>%
    unique()

write.table(mono.lps.2.eSNPs, "snps_lists/monocyte_lps_2hr.txt", quote=F, col.names=F, row.names=F, sep="\t")
write.table(mono.lps.24.eSNPs, "snps_lists/monocyte_lps_24hr.txt", quote=F, col.names=F, row.names=F, sep="\t")
write.table(mono.ifn.eSNPs, "snps_lists/monocyte_ifn.txt", quote=F, col.names=F, row.names=F, sep="\t")
write.table(mono.naive.eSNPs, "snps_lists/monocyte_naive.txt", quote=F, col.names=F, row.names=F, sep="\t")
