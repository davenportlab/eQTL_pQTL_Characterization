nextflow.enable.dsl = 2

params.genotypes = "/nfs/users/nfs_n/nm18/gains_team282/nikhil/data/genotypes/eigengene_sva_ss_genotypes.raw"
params.module_qtl_snps = "/nfs/users/nfs_n/nm18/gains_team282/nikhil/expression/eigengene_sva/mqtl_full_summary_statistics_snps.txt"
params.module_qtl_snps_all_pcs = "/nfs/users/nfs_n/nm18/gains_team282/nikhil/expression/eigengene_sva/mqtl_all_pcs_full_summary_statistics_snps.txt"
params.design_matrix = "/nfs/users/nfs_n/nm18/gains_team282/nikhil/expression/eigengene_sva/mapping_data.csv"
params.output_dir = "/nfs/users/nfs_n/nm18/gains_team282/nikhil/expression/eigengene_sva/wgcna_summary_statistics/"

process PERFORM_ASSOCIATION_TESTS {

    errorStrategy "retry"
    maxRetries 3

    label "R"

    publishDir "$params.output_dir/", mode: "move"

    input:
        val(me)

    output:
        path("*.tsv"), emit: module_associations
    
    script:

        """
        Rscript \\
            $workflow.projectDir/module_qtl_summary_statistics.R \\
            $params.genotypes \\
            $params.design_matrix \\
            $params.module_qtl_snps \\
            $me
        """
}

process PERFORM_ASSOCIATION_TESTS_ALL_PCS {

    errorStrategy "retry"
    maxRetries 3

    label "R"

    publishDir "$params.output_dir/", mode: "move"

    input:
        val(me)

    output:
        path("*.tsv"), emit: module_associations
    
    script:

        """
        Rscript \\
            $workflow.projectDir/module_qtl_summary_statistics_all_pcs.R \\
            $params.genotypes \\
            $params.design_matrix \\
            $params.module_qtl_snps_all_pcs \\
            $me
        """
}

workflow {

    me_channel = Channel
        .fromPath("$params.module_qtl_snps")
        .splitCsv(sep: "\t", skip: 1)
        .map{ row -> row[0] }
        .unique()
    
    me_all_pcs_channel = Channel
        .fromPath("$params.module_qtl_snps_all_pcs")
        .splitCsv(sep: "\t", skip: 1)
        .map{ row -> row[0] }
        .unique()

    PERFORM_ASSOCIATION_TESTS(me_channel)

    PERFORM_ASSOCIATION_TESTS_ALL_PCS(me_all_pcs_channel)
}
