nextflow.enable.dsl = 2

params.genotypes = "/nfs/users/nfs_n/nm18/gains_team282/nikhil/data/genotypes/eigengene_sva_genotypes.raw"
params.mapping_eigengenes = "/nfs/users/nfs_n/nm18/gains_team282/nikhil/expression/eigengene_sva/mapping_eigengenes.txt"
params.design_matrix = "/nfs/users/nfs_n/nm18/gains_team282/nikhil/expression/eigengene_sva/mapping_data.csv"
params.output_dir = "/nfs/users/nfs_n/nm18/gains_team282/nikhil/expression/eigengene_sva/initial_pass/"

process REGRESS_CIS_ESNPS {

    errorStrategy "retry"
    maxRetries 3

    label "R"

    output:
        path("regressed.eigengenes.RDS"), emit: regressed_eigengenes
    
    script:

        """
        Rscript \\
            $workflow.projectDir/eigengene_single_variant_association_regress_cis_esnps.R \\
            $params.genotypes \\
            $params.design_matrix
        """
}

process PERFORM_ASSOCIATION_TESTS {

    errorStrategy "retry"
    maxRetries 3

    label "R"

    publishDir "$params.output_dir/", mode: "move"

    input:
        val(me)
        path(regressed_eigengenes)

    output:
        path("${me}.tsv"), emit: module_associations
    
    script:

        """
        Rscript \\
            $workflow.projectDir/eigengene_single_variant_association.R \\
            $params.genotypes \\
            $params.design_matrix \\
            $regressed_eigengenes \\
            $me \\
            ${me}.tsv
        """
}

workflow {

    REGRESS_CIS_ESNPS()

    me_channel = Channel
        .fromPath("$params.mapping_eigengenes")
        .splitText()
        .map{ me -> me.strip() }

    PERFORM_ASSOCIATION_TESTS(
        me_channel,
        REGRESS_CIS_ESNPS.out.regressed_eigengenes
    )
}
