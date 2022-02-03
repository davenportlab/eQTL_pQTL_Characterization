args = commandArgs(trailingOnly = TRUE)
output.file = args[1]

cis.eqtl.summary <- readRDS("/lustre/scratch119/realdata/mdt3/projects/gains_team282/eqtl/cisresults/ciseqtl_all.rds")

write.table(unique(cis.eqtl.summary$snps), output.file, row.names=F, col.names=F, quote=F)
