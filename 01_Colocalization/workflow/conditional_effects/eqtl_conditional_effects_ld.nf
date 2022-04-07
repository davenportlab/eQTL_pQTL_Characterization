nextflow.enable.dsl = 2

params.chr = "1"
params.genotypes = "/nfs/users/nfs_n/nm18/gains_team282/Genotyping/All_genotyping_merged_filtered_b38_refiltered_rsID"
params.gene_exp = "/lustre/scratch119/humgen/projects/gains_team282/eqtl/data/logcpm_864_20412_hla.txt"
params.output_dir = "/nfs/users/nfs_n/nm18/gains_team282/nikhil/colocalization/cis_eqtl/conditional_effects/LD/"


process PREPARE_LOCI {

    label "Rbigmem"

    output:
        path("eqtl_individuals.txt"),       emit: mapping_individuals
        path("single_signal_stats.txt"),    emit: single_signal_stats
        path("ENSG*"),                      emit: locus_files

    script:

        """
        head -n 1 $params.gene_exp | sed 's/"//g' | sed 's/_.//g' | sed 's/[[:space:]]/\\n/g' | sed 's/[[:space:]]*\$//g' > gene_exp_individuals.txt
        cat ${params.genotypes}.fam | awk '{ print \$1; }' | grep -wFf gene_exp_individuals.txt | awk '{ print \$1 "\\t" \$1; }' > eqtl_individuals.txt

        plink \\
            --bfile $params.genotypes \\
            --keep eqtl_individuals.txt \\
            --freq \\
            --allow-extra-chr

        Rscript $workflow.projectDir/eqtl_conditional_effects_ld_prepare_loci.R $params.chr
        """
}

process RUN_COJO {

    errorStrategy "ignore"

    label "cojo"

    input:
        path(eqtl_individuals)
        tuple(val(locus), path("locus_files/*"))
    
    output:
        path("${locus}_conditional_signals.txt"), emit: conditional_signal_stats

    script:

        """
        for file_name in locus_files/*.cond.snps
        do

            gcta64 \\
                --bfile $params.genotypes \\
                --keep $eqtl_individuals \\
                --extract locus_files/${locus}.snps \\
                --cojo-file locus_files/${locus}.ma \\
                --cojo-cond \$file_name \\
                --cojo-collinear 0.99 \\
                --out \$(basename \$file_name)
        
        done

        for file_name in *.cma.cojo
        do

            gene_signal_pair=\$(echo \$file_name | sed 's/\\..*//g')
            gene=\$(echo \$gene_signal_pair | sed 's/-.*//g')
            signal=\$(echo \$gene_signal_pair | sed 's/.*-//g')
            awk -v gene=\$gene -v signal=\$signal 'BEGIN { OFS="\\t"; } NR > 1 { print gene, signal, \$1, \$2, \$3, \$11, \$12, \$13; }' \$file_name > \${file_name}.part
        
        done

        if compgen -G "*.part" > /dev/null; then
            cat *.part > ${locus}_conditional_signals.txt
        else
            touch ${locus}_conditional_signals.txt
        fi
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
        .map{ file -> [file.getSimpleName().replaceAll(/-.*/, ""), file] }
        .groupTuple()
    
    RUN_COJO(PREPARE_LOCI.out.mapping_individuals, locus_files)

    MERGE_CONDITIONAL_OUTPUT(
        PREPARE_LOCI.out.single_signal_stats,
        RUN_COJO.out.conditional_signal_stats.collect()
    )
}
