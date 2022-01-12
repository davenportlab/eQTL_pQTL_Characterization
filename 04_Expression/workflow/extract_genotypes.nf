nextflow.enable.dsl = 2

params.genotypes_prefix = "~/gains_team282/Genotyping/All_genotyping_merged_filtered_b38_refiltered_rsID"
params.mapping_patients = "~/gains_team282/nikhil/expression/eigengene_sva/mapping_patients.txt"
params.output_dir = "/nfs/users/nfs_n/nm18/gains_team282/nikhil/expression/eigengene_sva/genotypes/"

process EXTRACT_GENOTYPES {

    errorStrategy "retry"
    maxRetries 3

    label "plink"

    publishDir "$params.output_dir/", mode: "move"

    output:
        path("*.raw"), emit: genotype_data_for_mapping

    script:

        """
        for CHR in {1..23}
        do
            plink \\
                --bfile $params.genotypes_prefix \\
                --keep $params.mapping_patients \\
                --recode A \\
                --out ./eigengene_sva_genotypes_\${CHR} \\
                --allow-extra-chr \\
                --chr \${CHR} \\
                --maf 0.01
        done
        """
}

workflow {

    EXTRACT_GENOTYPES()
}
