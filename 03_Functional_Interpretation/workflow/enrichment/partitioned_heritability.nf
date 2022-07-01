nextflow.enable.dsl = 2

params.grm_dir = "/nfs/users/nfs_n/nm18/gains_team282/epigenetics/enrichment/heritability/"
params.eigengenes = "/nfs/users/nfs_n/nm18/gains_team282/nikhil/expression/eigengene_sva/mapping_eigengenes.txt"
params.output_dir = "/nfs/users/nfs_n/nm18/gains_team282/epigenetics/enrichment/partitioned_heritability/"


process ESTIMATE_PARTITIONED_HERITABILITY {

    label "heritability"

    input:
        val(annotation)

    output:
        path("*.variance_components.csv"), emit: variance_components

    script:
    
        """
        Rscript $workflow.projectDir/partitioned_heritability.R $annotation
        """
}

process COLLATE_RESULTS {

    label "simple_bash"

    publishDir "$params.output_dir/", mode: "move"

    input:
        path("var_comps/*")
    
    output:
        path("variance_components.csv")
    
    script:

        """
        echo -e "Eigengene\\tAnnotation\\tComponent\\tVariance\\tProportion" > variance_components.csv

        cat var_comps/*.csv >> variance_components.csv
        """
}

workflow {

    grm_channel = Channel
        .fromPath("$params.grm_dir/*", type: "dir")
        .map{ item -> item.getSimpleName() }

    ESTIMATE_PARTITIONED_HERITABILITY(grm_channel)

    COLLATE_RESULTS(ESTIMATE_PARTITIONED_HERITABILITY.out.variance_components.collect())
}
