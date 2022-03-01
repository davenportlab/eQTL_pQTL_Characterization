nextflow.enable.dsl = 2

params.lead_snps = "/nfs/users/nfs_n/nm18/gains_team282/nikhil/colocalization/cis_eqtl/fine_mapping/LD/lead_snps.80r2.tags.tsv"
params.conditional_snps = "/nfs/users/nfs_n/nm18/gains_team282/nikhil/colocalization/cis_eqtl/fine_mapping/LD/conditional_snps.80r2.tags.tsv"
params.output_dir = "/nfs/users/nfs_n/nm18/gains_team282/epigenetics/calderon_et_al_hg19/"


process PREPARE_SNP_LIST {

    label "simple_bash"

    output:
        path("lead_snps.txt"),          emit: lead_snps
        path("conditional_snps.txt"),   emit: conditional_snps
    
    
}


workflow {


}
