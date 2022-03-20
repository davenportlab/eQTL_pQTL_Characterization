args = commandArgs(trailingOnly = TRUE)
cis.eqtl.snps.output.file = args[1]
cis.pqtl.snps.output.file = args[2]
trans.pqtl.snps.output.file = args[3]
gene.exp.sample.file = args[4]
protein.exp.sample.file = args[5]

cis.eqtl.summary <- readRDS("/lustre/scratch119/realdata/mdt3/projects/gains_team282/eqtl/cisresults/ciseqtl_all.rds")

write.table(unique(cis.eqtl.summary$snps), cis.eqtl.snps.output.file, row.names=F, col.names=F, quote=F)

cis.pqtl <- readRDS("/nfs/users/nfs_n/nm18/gains_team282/proteomics/pqtl/pqtl_ms2019/cis_pqtl_all.RDS")
trans.pqtl <- readRDS("/nfs/users/nfs_n/nm18/gains_team282/proteomics/pqtl/pqtl_ms2019/trans_pqtl_all.RDS")

write.table(unique(cis.pqtl$SNP), cis.pqtl.snps.output.file, row.names=F, col.names=F, quote=F)
write.table(unique(trans.pqtl$SNP), trans.pqtl.snps.output.file, row.names=F, col.names=F, quote=F)

# Load gene/protein expression to get sample names
gene.exp <- read.table("/lustre/scratch119/humgen/projects/gains_team282/eqtl/data/logcpm_864_20412_hla.txt")
protein.sample.info <- read.csv("/nfs/users/nfs_n/nm18/gains_team282/proteomics/MS2019_processed_data/sample_info_1860_MS2019.csv")

gene.exp.samples <- unique(gsub("_.*$", "", colnames(gene.exp)))
gene.exp.samples <- cbind(gene.exp.samples, gene.exp.samples)
write.table(gene.exp.samples, gene.exp.sample.file, row.names=F, col.names=F, quote=F, sep="\t")

protein.exp.samples <- unique(protein.sample.info$Patient)
protein.exp.samples <- cbind(protein.exp.samples, protein.exp.samples)
write.table(protein.exp.samples, protein.exp.sample.file, row.names=F, col.names=F, quote=F, sep="\t")
