nextflow.enable.dsl = 2

params.output_dir = "/nfs/users/nfs_n/nm18/gains_team282/nikhil/colocalization/cis_eqtl/fine_mapping/"
params.chr = "1"

process SPLIT_LOCI {

    label "Rbigmem"

    output:
        path("*.summary.RDS"),      emit: loci
        path("*.genotypes.RDS"),    emit: genotypes
        path("*.var_y.RDS"),        emit: var_y
        path("n_samples.RDS"),      emit: n_samples
    
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
        path(genotypes)
        path(var_y)
        path(n_samples)

    output:
        path("*.tsv"), emit: credible_set

    script:

        """
        Rscript $workflow.projectDir/susie_fine_mapping.R ${locus.getSimpleName()} $locus $genotypes $var_y $n_samples ${locus.getSimpleName()}.tsv
        """
}

process AGGREGATE_CREDIBILE_SETS {

    errorStrategy "retry"
    maxRetries 3

    publishDir "$params.output_dir/", mode: "move"

    input:
        path("credible_sets/*.tsv")

    output:
        path("chr${params.chr}_credible_sets.tsv")

    script:

        """
        echo -e "Gene\tSNP\tSNP_Prob\tCredible_Set" > chr${params.chr}_credible_sets.tsv

        cat credible_sets/*.tsv >> chr${params.chr}_credible_sets.tsv
        """
}


workflow {

    SPLIT_LOCI()

    SUSIE_FINE_MAPPING(
        SPLIT_LOCI.out.loci.flatten(),
        SPLIT_LOCI.out.genotypes.flatten(),
        SPLIT_LOCI.out.var_y.flatten(),
        SPLIT_LOCI.out.n_samples
    )

    AGGREGATE_CREDIBILE_SETS(SUSIE_FINE_MAPPING.out.credible_set.collect())
}
