nextflow.enable.dsl = 2

params.genotypes_prefix = "/nfs/users/nfs_n/nm18/gains_team282/Genotyping/All_genotyping_merged_filtered_b38_refiltered_rsID"
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
        Rscript $workflow.projectDir/extract_genotypes.R cis.eqtl.snps.txt cis.pqtl.snps.txt trans.pqtl.snps.txt gene.exp.samples.txt protein.exp.samples.txt

        plink \\
            --bfile $params.genotypes_prefix \\
            --keep protein.exp.samples.txt \\
            --extract cis.pqtl.snps.txt \\
            --recode A \\
            --out ./cis_pqtl_genotypes \\
            --allow-extra-chr

        plink \\
            --bfile $params.genotypes_prefix \\
            --keep protein.exp.samples.txt \\
            --extract trans.pqtl.snps.txt \\
            --recode A \\
            --out ./trans_pqtl_genotypes \\
            --allow-extra-chr

        for CHR in {1..22}
        do  
            plink \\
                --bfile $params.genotypes_prefix \\
                --keep gene.exp.samples.txt \\
                --extract cis.eqtl.snps.txt \\
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
