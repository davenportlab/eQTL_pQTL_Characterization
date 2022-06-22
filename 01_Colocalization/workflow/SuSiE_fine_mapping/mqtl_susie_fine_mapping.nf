nextflow.enable.dsl = 2

params.output_dir = "/nfs/users/nfs_n/nm18/gains_team282/nikhil/colocalization/mqtl/fine_mapping/SuSiE/"

process SPLIT_LOCI {

    label "Rbigmem"

    output:
        path("*.summary.RDS"),      emit: susie_loci
        path("*.genotypes.RDS"),    emit: susie_genotypes
        path("*.var.y.RDS"),        emit: susie_var_y
        path("n_samples.RDS"),      emit: n_samples
    
    script:

        """
        Rscript $workflow.projectDir/mqtl_susie_fine_mapping_split_loci.R
        """
}

process SUSIE_FINE_MAPPING {

    errorStrategy "retry"
    maxRetries 1

    label "R_cor"

    input:
        path(locus)
        path(genotypes)
        path(var_y)
        path(n_samples)

    output:
        path("${locus.getSimpleName()}.tsv"), emit: credible_set

    script:

        """
        Rscript $workflow.projectDir/mqtl_susie_fine_mapping.R \\
            ${locus.getSimpleName()} \\
            $locus \\
            $genotypes \\
            $var_y \\
            $n_samples \\
            10 \\
            ${locus.getSimpleName()}.tsv
        """
}

process AGGREGATE_CREDIBILE_SETS {

    errorStrategy "retry"
    maxRetries 3

    publishDir "$params.output_dir/", mode: "move"

    input:
        path("cs/*.tsv")

    output:
        path("module_qtl_credible_sets.tsv")

    script:

        """
        echo -e "Gene\tSNP\tSNP_Prob\tCredible_Set\tNotes" > module_qtl_credible_sets.tsv
        cat cs/*.tsv >> module_qtl_credible_sets.tsv
        """
}


workflow {

    SPLIT_LOCI()

    SUSIE_FINE_MAPPING(
        SPLIT_LOCI.out.susie_loci.flatten(),
        SPLIT_LOCI.out.susie_genotypes.flatten(),
        SPLIT_LOCI.out.susie_var_y.flatten(),
        SPLIT_LOCI.out.n_samples
    )

    AGGREGATE_CREDIBILE_SETS(
        SUSIE_FINE_MAPPING.out.credible_set.collect()
    )
}
