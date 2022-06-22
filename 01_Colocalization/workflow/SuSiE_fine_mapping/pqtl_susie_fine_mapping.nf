nextflow.enable.dsl = 2

params.output_dir = "/nfs/users/nfs_n/nm18/gains_team282/nikhil/colocalization/pqtl/fine_mapping/SuSiE/"

process SPLIT_LOCI {

    label "Rbigmem"

    output:
        path("*.summary.RDS"),      emit: susie_loci
        path("*.genotypes.RDS"),    emit: susie_genotypes
        path("*.var.y.RDS"),        emit: susie_var_y
        path("n_samples.RDS"),      emit: n_samples
    
    script:

        """
        Rscript $workflow.projectDir/pqtl_susie_fine_mapping_split_loci.R
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
        path("${locus.getSimpleName()}.tsv"), emit: credible_set

    script:

        """
        Rscript $workflow.projectDir/pqtl_susie_fine_mapping.R \\
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
        path("cis_cs/*.tsv")
        path("trans_cs/*.tsv")

    output:
        path("cis_pqtl_credible_sets.tsv")
        path("trans_pqtl_credible_sets.tsv")

    script:

        """
        echo -e "Gene\tSNP\tSNP_Prob\tCredible_Set\tNotes" > cis_pqtl_credible_sets.tsv
        cat cis_cs/*.tsv >> cis_pqtl_credible_sets.tsv

        echo -e "pQTL_ID\tSNP\tSNP_Prob\tCredible_Set\tNotes" > trans_pqtl_credible_sets.tsv
        cat trans_cs/*.tsv >> trans_pqtl_credible_sets.tsv
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

    cs_channel = SUSIE_FINE_MAPPING.out.credible_set
        .branch {
            cis: it =~ /ENSG/
            trans: !(it =~ /ENSG/)
        }

    AGGREGATE_CREDIBILE_SETS(
        cs_channel.cis.collect(),
        cs_channel.trans.collect()
    )
}
