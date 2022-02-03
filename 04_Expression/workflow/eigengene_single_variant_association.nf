nextflow.enable.dsl = 2

params.genotypes_dir = "/nfs/users/nfs_n/nm18/gains_team282/nikhil/data/genotypes/"
params.design_matrix = "/nfs/users/nfs_n/nm18/gains_team282/nikhil/expression/eigengene_sva/mapping_data.csv"
params.output_dir = "/nfs/users/nfs_n/nm18/gains_team282/nikhil/expression/eigengene_sva/initial_pass/"
params.me = "ME_1"

process PERFORM_ASSOCIATION_TESTS {

    errorStrategy "retry"
    maxRetries 3

    label "R"

    input:
        path(genotypes_file)

    output:
        path("*.tsv"), emit: module_associations
    
    script:

        """
        CHR=\$(basename $genotypes_file .raw | cut -c25-26)

        Rscript \\
            $workflow.projectDir/eigengene_single_variant_association.R \\
            $genotypes_file \\
            $params.design_matrix \\
            $params.me \\
            ${params.me}_\${CHR}.tsv
        """
}

process AGGREGATE_ASSOCIATION_RESULTS {

    errorStrategy "retry"
    maxRetries 3

    label "multi_cpu_bash"

    publishDir "$params.output_dir/", mode: "move"

    input:
        path("associations/*.tsv")
    
    output:
        path("${params.me}.tsv")

    script:

        """
        cat associations/*.tsv > ${params.me}.tsv
        """
}

workflow {

    geno_channel = Channel.fromPath("${params.genotypes_dir}/*.raw")

    PERFORM_ASSOCIATION_TESTS(geno_channel)

    AGGREGATE_ASSOCIATION_RESULTS(PERFORM_ASSOCIATION_TESTS.out.module_associations.collect())
}
