nextflow.enable.dsl = 2

params.genotypes_prefix = "/nfs/users/nfs_n/nm18/gains_team282/Genotyping/All_genotyping_merged_filtered_b38_refiltered_rsID"
params.mapping_patients = "/nfs/users/nfs_n/nm18/gains_team282/nikhil/expression/eigengene_sva/mapping_patients.txt"
params.mqtl_snps = "/nfs/users/nfs_n/nm18/gains_team282/nikhil/expression/eigengene_sva/mqtl_full_summary_statistics_snps.txt"
params.mqtl_pc_snps = "/nfs/users/nfs_n/nm18/gains_team282/nikhil/expression/eigengene_sva/mqtl_all_pcs_full_summary_statistics_snps.txt"
params.replication_snps = "/nfs/users/nfs_n/nm18/gains_team282/nikhil/expression/eigengene_sva/all_mqtl_all_pcs.csv"
params.output_dir = "/nfs/users/nfs_n/nm18/gains_team282/nikhil/data/genotypes/"

process EXTRACT_GENOTYPES {

    errorStrategy "retry"
    maxRetries 3

    label "plink"

    publishDir "$params.output_dir/", mode: "move"

    output:
        path("*.raw")

    script:

        """
        awk 'NR > 1 { print \$2; }' $params.mqtl_snps > mqtl.snp_set.txt
        awk 'NR > 1 { print \$2; }' $params.mqtl_pc_snps > mqtl.snp_set_all_pcs.txt
        cat mqtl.snp_set.txt mqtl.snp_set_all_pcs.txt | sort | uniq > mqtl.snps.txt

        plink \\
            --bfile $params.genotypes_prefix \\
            --keep $params.mapping_patients \\
            --extract mqtl.snps.txt \\
            --recode A \\
            --out ./eigengene_sva_ss_genotypes \\
            --allow-extra-chr \\
            --maf 0.01
        """
}

process EXTRACT_REPLICATION_GENOTYPES() {

    errorStrategy "retry"
    maxRetries 3

    label "plink"

    publishDir "$params.output_dir/", mode: "move"

    output:
        path("*.raw")

    script:

        """
        awk -F ',' 'NR > 1 { print \$1; }' $params.replication_snps | sort | uniq > rep.snps.txt

        plink \\
            --bfile $params.genotypes_prefix \\
            --extract rep.snps.txt \\
            --recode A \\
            --out ./eigengene_sva_rep_genotypes \\
            --allow-extra-chr \\
            --maf 0.01
        """
}

workflow {

    EXTRACT_GENOTYPES()

    EXTRACT_REPLICATION_GENOTYPES()
}
