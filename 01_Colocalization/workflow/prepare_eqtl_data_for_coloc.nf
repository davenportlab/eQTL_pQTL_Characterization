nextflow.enable.dsl = 2

params.vcf_dir = "/lustre/scratch118/humgen/resources/1000g/release/20201028"
params.eur_samples = "~/eQTL_pQTL_Characterization/01_Colocalization/data/1000G/EUR.samples.txt"
params.output_dir = "~/gains_team282/nikhil/colocalization/cis_eqtl/"

process GENERATE_LOCUS_DATA {

    errorStrategy "retry"
    maxRetries 5

    label "R"

    output:
        path("*.csv"), emit: locus_data_for_ld_calculation
        path("*.csv"), emit: locus_data_for_aggregation

    script:

        """
        Rscript $workflow.projectDir/prepare_eqtl_data_for_coloc_split.R
        """
}

process CALCULATE_LD {

    errorStrategy "retry"
    maxRetries 5

    label "vcftools"

    input:
        path(locus_file)

    output:
        path("*_ld.csv"), emit: locus_ld

    script:

        """
        CHR=\$(awk -F ',' 'NR==2 { print \$5; }' $locus_file)
        VCF_FILE=\$(ls /lustre/scratch118/humgen/resources/1000g/release/20201028 | grep "chr\${CHR}.*.gz\$")
        NAME=\$(basename $locus_file .csv)

        awk -F ',' 'NR!=1 { print "chr" \$5 "\t" \$2; }' $locus_file > snps.txt

        bcftools view \\
            -R snps.txt \\
            $params.vcf_dir/\$VCF_FILE | \\
            $workflow.projectDir/prepare_eqtl_data_for_coloc_ld_matrix/filter_snps.py $locus_file | \\
            vcftools \\
            --vcf - \\
            --keep $params.eur_samples \\
            --mac 1 \\
            --geno-r2 \\
            --stdout | \\
            $workflow.projectDir/prepare_eqtl_data_for_coloc_ld_matrix/ld_matrix snps.txt > \${NAME}_ld.csv
        """
}

process AGGREGATE {

    errorStrategy "retry"
    maxRetries 5

    label "R"

    input:
        path("locus_data/*")
        path("locus_ld/*")

    script:

        """
        Rscript $workflow.projectDir/prepare_eqtl_data_for_coloc_aggregate.R
        """
}


workflow {

    GENERATE_LOCUS_DATA()

    CALCULATE_LD(GENERATE_LOCUS_DATA.out.locus_data_for_ld_calculation.buffer(size: 20))

    AGGREGATE(
        GENERATE_LOCUS_DATA.out.locus_data_for_aggregation,
        CALCULATE_LD.out.locus_ld.collect().flatten()
    )
}
