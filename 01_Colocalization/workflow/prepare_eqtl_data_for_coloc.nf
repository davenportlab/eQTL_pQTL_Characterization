nextflow.enable.dsl = 2

params.vcf_dir = "/lustre/scratch118/humgen/resources/1000g/release/20201028"
params.eur_samples = "/nfs/users/nfs_n/nm18/eQTL_pQTL_Characterization/01_Colocalization/data/1000G/EUR.samples.txt"
params.output_dir = "/nfs/users/nfs_n/nm18/gains_team282/nikhil/colocalization/cis_eqtl/"
params.chr = "1"

process GENERATE_LOCUS_DATA {

    errorStrategy "retry"
    maxRetries 5

    label "Rbigmem"

    output:
        path("*.csv"), emit: locus_data_for_ld_calculation
        path("*.csv"), emit: locus_data_for_aggregation

    script:

        """
        Rscript $workflow.projectDir/prepare_eqtl_data_for_coloc_split.R $params.chr
        """
}

process CALCULATE_LD {

    errorStrategy "retry"
    maxRetries 5

    label "vcftools"

    input:
        path(locus_file)

    output:
        path("*_ld.tsv"), emit: locus_ld

    script:

        """
        VCF_FILE=\$(ls $params.vcf_dir | grep "chr${params.chr}.filtered.*.gz\$")
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
            --stdout > ld.tsv || true

        awk 'NR!=1 { print \$2; }' ld.tsv > snp_col_1.txt
        awk 'NR!=1 { print \$3; }' ld.tsv > snp_col_2.txt
        cat snp_col_1.txt snp_col_2.txt | sort -u > snps.filtered.txt

        cat ld.tsv | $workflow.projectDir/prepare_eqtl_data_for_coloc_ld_matrix/ld_matrix snps.filtered.txt > \${NAME}_ld.tsv

        rm ld.tsv
        rm snp_col_1.txt
        rm snp_col_2.txt
        rm snps.txt
        rm snps.filtered.txt
        """
}

process AGGREGATE {

    errorStrategy "retry"
    maxRetries 5

    label "Rbigmem"

    input:
        path("locus_data/*")
        path("locus_ld/*")

    script:

        """
        Rscript $workflow.projectDir/prepare_eqtl_data_for_coloc_aggregate.R $params.chr
        """
}


workflow {

    GENERATE_LOCUS_DATA()

    CALCULATE_LD(GENERATE_LOCUS_DATA.out.locus_data_for_ld_calculation.flatten())

    AGGREGATE(
        GENERATE_LOCUS_DATA.out.locus_data_for_aggregation,
        CALCULATE_LD.out.locus_ld.collect()
    )
}
