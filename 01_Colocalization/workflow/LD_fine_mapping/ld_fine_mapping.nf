nextflow.enable.dsl = 2

params.genotypes_prefix = "/nfs/users/nfs_n/nm18/gains_team282/Genotyping/All_genotyping_merged_filtered_b38_refiltered_rsID"
params.output_dir = "/nfs/users/nfs_n/nm18/gains_team282/nikhil/colocalization/cis_eqtl/fine_mapping/LD/"
params.chr = "1"


process RETRIEVE_SNPS {

    label "Rbigmem"

    errorStrategy "retry"
    maxRetries 3

    output:
        path("*.txt"), emit: snps_files

    script:

        """
        Rscript $workflow.projectDir/ld_fine_mapping_retrieve_snps.R $params.chr
        """
}

process LD_FINE_MAPPING {

    label "plink"

    errorStrategy "retry"
    maxRetries 3

    publishDir "$params.output_dir/", mode: "move"

    input:
        path(snps_file)

    output:
        path("*.tags.tsv")
    
    script:

        """
        plink \\
            --bfile $params.genotypes_prefix \\
            --show-tags $snps_file \\
            --allow-extra-chr \\
            --tag-kb 1000 \\
            --tag-r2 1 \\
            --list-all
        
        cat plink.tags.list | sed -e "s/[[:space:]]\\+/\\t/g" | sed -e "s/^\\t//g" > ${snps_file.getSimpleName()}.100r2.tags.tsv

        plink \\
            --bfile $params.genotypes_prefix \\
            --show-tags $snps_file \\
            --allow-extra-chr \\
            --tag-kb 1000 \\
            --tag-r2 0.9 \\
            --list-all
        
        cat plink.tags.list | sed -e "s/[[:space:]]\\+/\\t/g" | sed -e "s/^\\t//g" > ${snps_file.getSimpleName()}.90r2.tags.tsv

        plink \\
            --bfile $params.genotypes_prefix \\
            --show-tags $snps_file \\
            --allow-extra-chr \\
            --tag-kb 1000 \\
            --tag-r2 0.8 \\
            --list-all
        
        cat plink.tags.list | sed -e "s/[[:space:]]\\+/\\t/g" | sed -e "s/^\\t//g" > ${snps_file.getSimpleName()}.80r2.tags.tsv
        """
}


workflow {

    RETRIEVE_SNPS()

    LD_FINE_MAPPING(RETRIEVE_SNPS.out.snps_files.flatten())
}
