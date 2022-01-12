nextflow.enable.dsl = 2

params.output_dir = "~/gains_team282/nikhil/colocalization/cis_eqtl/"
params.chr = "1"

process SPLIT_LOCI {

    label "Rbigmem"

    output:
        path("*.RDS"), emit: loci
    
    script:

        """
        Rscript $workflow.projectDir/susie_fine_mapping_split_loci.R $params.chr
        """
}

process SUSIE_FINE_MAPPING {

    errorStrategy "retry"
    maxRetries 3

    label "R"

    input:
        path(locus)

    output:
        path("*.tsv"), emit: credible_set

    script:

        """
        Rscript $workflow.projectDir/susie_fine_mapping.R ${locus.getSimpleName()} $locus ${locus.getSimpleName()}.tsv
        """
}

process AGGREGATE_CREDIBILE_SETS {

    errorStrategy "retry"
    maxRetries 3

    publishDir "$params.output_dir/", mode: "move"

    input:
        path("credible_sets/*.tsv")

    output:
        path("credible_sets.tsv")

    script:

        """
        echo -e "Gene\tSNP\tSNP_Prob\tCredible_Set" > credible_sets.tsv

        cat credible_sets/*.tsv >> credible_sets.tsv
        """
}


workflow {

    SPLIT_LOCI()

    SUSIE_FINE_MAPPING(SPLIT_LOCI.out.loci.flatten())

    AGGREGATE_CREDIBILE_SETS(SUSIE_FINE_MAPPING.out.credible_set.collect())
}
