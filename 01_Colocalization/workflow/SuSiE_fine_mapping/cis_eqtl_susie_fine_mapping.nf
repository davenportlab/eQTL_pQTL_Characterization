nextflow.enable.dsl = 2

params.output_dir = "/nfs/users/nfs_n/nm18/gains_team282/nikhil/colocalization/cis_eqtl/fine_mapping/SuSiE/"
params.chr = "1"

process SPLIT_LOCI {

    label "Rbigmem"

    output:
        path("full/*.summary.RDS"),             emit: full_susie_loci
        path("full/*.genotypes.RDS"),           emit: full_susie_genotypes
        path("full/*.var.y.RDS"),               emit: full_susie_var_y
        path("n_samples.RDS"),                  emit: full_n_samples
        path("conditional/*.summary.RDS"),      emit: conditional_susie_loci
        path("conditional/*.genotypes.RDS"),    emit: conditional_susie_genotypes
        path("conditional/*.var.y.RDS"),        emit: conditional_susie_var_y
        path("n_samples.RDS"),                  emit: conditional_n_samples
    
    script:

        """
        mkdir full/
        mkdir conditional/

        Rscript $workflow.projectDir/cis_eqtl_susie_fine_mapping_split_loci.R $params.chr
        """
}

process SUSIE_FINE_MAPPING_FULL {

    errorStrategy "retry"
    maxRetries 3

    label "R"

    input:
        path(locus)
        path(genotypes)
        path(var_y)
        path(n_samples)

    output:
        path("${locus.getSimpleName()}.tsv"), emit: credible_set

    script:

        """
        Rscript $workflow.projectDir/cis_eqtl_susie_fine_mapping.R \\
            ${locus.getSimpleName()} \\
            $locus \\
            $genotypes \\
            $var_y \\
            $n_samples \\
            10 \\
            ${locus.getSimpleName()}.tsv
        """
}

process SUSIE_FINE_MAPPING_CONDITIONAL {

    errorStrategy "retry"
    maxRetries 3

    label "R"

    input:
        path(locus)
        path(genotypes)
        path(var_y)
        path(n_samples)

    output:
        path("${locus.getSimpleName()}.tsv"), emit: credible_set

    script:

        """
        Rscript $workflow.projectDir/cis_eqtl_susie_fine_mapping.R \\
            ${locus.getSimpleName()} \\
            $locus \\
            $genotypes \\
            $var_y \\
            $n_samples \\
            1 \\
            ${locus.getSimpleName()}.tsv
        """
}

process AGGREGATE_CREDIBILE_SETS {

    errorStrategy "retry"
    maxRetries 3

    publishDir "$params.output_dir/", mode: "move"

    input:
        path("full_cs/*.tsv")
        path("conditional_cs/*.tsv")

    output:
        path("full/full_chr${params.chr}_credible_sets.tsv")
        path("conditional/conditional_chr${params.chr}_credible_sets.tsv")

    script:

        """
        mkdir full/
        mkdir conditional/

        echo -e "Gene\tSNP\tSNP_Prob\tCredible_Set\tNotes" > full/full_chr${params.chr}_credible_sets.tsv
        cat full_cs/*.tsv >> full/full_chr${params.chr}_credible_sets.tsv

        echo -e "Gene\tSNP\tSNP_Prob\tCredible_Set\tNotes" > conditional/conditional_chr${params.chr}_credible_sets.tsv
        cat conditional_cs/*.tsv >> conditional/conditional_chr${params.chr}_credible_sets.tsv
        """
}


workflow {

    SPLIT_LOCI()

    SUSIE_FINE_MAPPING_FULL(
        SPLIT_LOCI.out.full_susie_loci.flatten(),
        SPLIT_LOCI.out.full_susie_genotypes.flatten(),
        SPLIT_LOCI.out.full_susie_var_y.flatten(),
        SPLIT_LOCI.out.full_n_samples
    )

    SUSIE_FINE_MAPPING_CONDITIONAL(
        SPLIT_LOCI.out.conditional_susie_loci.flatten(),
        SPLIT_LOCI.out.conditional_susie_genotypes.flatten(),
        SPLIT_LOCI.out.conditional_susie_var_y.flatten(),
        SPLIT_LOCI.out.conditional_n_samples
    )

    AGGREGATE_CREDIBILE_SETS(
        SUSIE_FINE_MAPPING_FULL.out.credible_set.collect(),
        SUSIE_FINE_MAPPING_CONDITIONAL.out.credible_set.collect()
    )
}
