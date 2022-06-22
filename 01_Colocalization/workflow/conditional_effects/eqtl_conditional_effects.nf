nextflow.enable.dsl = 2

params.chr = "1"
params.genotypes = "/nfs/users/nfs_n/nm18/gains_team282/nikhil/data/genotypes/eqtl_genotypes_1.raw"
params.output_dir = "/nfs/users/nfs_n/nm18/gains_team282/nikhil/colocalization/cis_eqtl/conditional_effects/LMM/"


process PREPARE_LOCI {

    label "Rbigmem"

    output:
        path("single_signal_stats.txt"),    emit: single_signal_stats
        path("design_matrix.csv"),          emit: design_matrix
        path("ENSG*"),                      emit: locus_files

    script:

        """
        Rscript $workflow.projectDir/eqtl_conditional_effects_prepare_loci.R $params.chr
        """
}

process RUN_LMM {

    errorStrategy "retry"
    maxRetries 3

    label "lmm"

    input:
        path(design_matrix)
        tuple(val(locus), path("locus_files/*"))
    
    output:
        path("${locus}-*.tsv"), emit: conditional_signal_stats

    script:

        """
        Rscript $workflow.projectDir/eqtl_conditional_effects_association.R \\
            $params.genotypes \\
            $design_matrix \\
            $locus \\
            locus_files/${locus}.snps.txt \\
            locus_files/${locus}.conditional_snps.txt
        """
}

process MERGE_CONDITIONAL_OUTPUT {

    label "simple_bash"

    publishDir "$params.output_dir/", mode: "move"

    input:
        path(single_signal_stats)
        path("conditional_signals/*")
    
    output:
        path("chr${params.chr}_conditional_cis_eQTL_summary_statistics.tsv")
    
    script:

        """
        echo -e "Gene\tSignal\tChr\tSNP\tPosition\tBeta\tSE\tP_Value" > chr${params.chr}_conditional_cis_eQTL_summary_statistics.tsv

        cat $single_signal_stats conditional_signals/* >> chr${params.chr}_conditional_cis_eQTL_summary_statistics.tsv
        """
}


workflow {

    PREPARE_LOCI()

    locus_files = PREPARE_LOCI.out.locus_files
        .flatten()
        .map{ file -> [file.getSimpleName().replaceAll(/\..*/, ""), file] }
        .groupTuple()
    
    RUN_LMM(PREPARE_LOCI.out.design_matrix, locus_files)

    MERGE_CONDITIONAL_OUTPUT(
        PREPARE_LOCI.out.single_signal_stats,
        RUN_LMM.out.conditional_signal_stats.collect()
    )
}
