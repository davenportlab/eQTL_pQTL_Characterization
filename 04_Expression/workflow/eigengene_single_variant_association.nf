nextflow.enable.dsl = 2

params.genotypes = "/nfs/users/nfs_n/nm18/gains_team282/nikhil/data/genotypes/eigengene_sva_genotypes.raw"
params.mapping_eigengenes = "/nfs/users/nfs_n/nm18/gains_team282/nikhil/expression/eigengene_sva/mapping_eigengenes.txt"
params.design_matrix = "/nfs/users/nfs_n/nm18/gains_team282/nikhil/expression/eigengene_sva/mapping_data.csv"
params.output_dir = "/nfs/users/nfs_n/nm18/gains_team282/nikhil/expression/eigengene_sva/initial_pass/"

process PERFORM_ASSOCIATION_TESTS {

    errorStrategy "retry"
    maxRetries 3

    label "R"

    publishDir "$params.output_dir/", mode: "move"

    input:
        val(me)

    output:
        path("${me}.tsv"), emit: module_associations
    
    script:

        """
        Rscript \\
            $workflow.projectDir/eigengene_single_variant_association.R \\
            $params.genotypes \\
            $params.design_matrix \\
            $me \\
            ${me}.tsv
        """
}

workflow {

    me_channel = Channel
        .fromPath("$params.mapping_eigengenes")
        .splitText()
        .map{ me -> me.strip() }

    PERFORM_ASSOCIATION_TESTS(me_channel)
}
