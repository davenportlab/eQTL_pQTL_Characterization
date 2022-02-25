nextflow.enable.dsl = 2

params.genotypes_prefix = "/nfs/users/nfs_n/nm18/gains_team282/Genotyping/All_genotyping_merged_filtered_b38_refiltered_rsID"
params.mapping_patients = "/nfs/users/nfs_n/nm18/gains_team282/nikhil/expression/eigengene_sva/mapping_patients.txt"
params.mqtl_snps = "/nfs/users/nfs_n/nm18/gains_team282/nikhil/expression/eigengene_sva/mqtl_snp_table.csv"
params.output_dir = "/nfs/users/nfs_n/nm18/gains_team282/nikhil/data/genotypes/"

process EXTRACT_GENOTYPES {

    errorStrategy "retry"
    maxRetries 3

    label "plink"

    publishDir "$params.output_dir/", mode: "move"

    output:
        path("*.raw"), emit: genotype_data_for_mapping

    script:

        """
        Rscript $workflow.projectDir/extract_genotypes.R eqtl.snps.txt

        awk -F "," 'NR > 1 { print \$1; }' $params.mqtl_snps > mqtl.snps.txt

        plink \\
            --bfile $params.genotypes_prefix \\
            --keep $params.mapping_patients \\
            --extract mqtl.snps.txt \\
            --recode A \\
            --out ./eigengene_sva_genotypes \\
            --allow-extra-chr \\
            --maf 0.01

        for CHR in {1..22}
        do  
            plink \\
                --bfile $params.genotypes_prefix \\
                --keep $params.mapping_patients \\
                --extract eqtl.snps.txt \\
                --recode A \\
                --out ./eqtl_genotypes_\${CHR} \\
                --allow-extra-chr \\
                --chr \${CHR}
        done
        """
}

workflow {

    EXTRACT_GENOTYPES()
}
