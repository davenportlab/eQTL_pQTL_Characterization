nextflow.enable.dsl = 2

params.genotypes = "/nfs/users/nfs_n/nm18/gains_team282/Genotyping/All_genotyping_merged_filtered_b38_refiltered_rsID"
params.protein_sample_info = "/nfs/users/nfs_n/nm18/gains_team282/proteomics/MS2019_processed_data/sample_info_1860_MS2019.csv"
params.output_dir = "/nfs/users/nfs_n/nm18/gains_team282/nikhil/colocalization/pqtl/conditional_effects/FINEMAP/"


process PREPARE_LOCI {

    label "Rbigmem"

    output:
        path("pqtl_individuals.txt"),           emit: mapping_individuals
        path("cis_single_signal_stats.txt"),    emit: cis_single_signal_stats
        path("ENSG*"),                          emit: cis_locus_files
        path("trans_single_signal_stats.txt"),  emit: trans_single_signal_stats
        path("Trans*"),                         emit: trans_locus_files

    script:

        """
        sed 's/"//g' $params.protein_sample_info | awk -F "," 'NR > 1 { print \$4; }' > protein_exp_individuals.txt
        cat ${params.genotypes}.fam | awk '{ print \$1; }' | grep -wFf protein_exp_individuals.txt | awk '{ print \$1 "\\t" \$1; }' > pqtl_individuals.txt

        plink \\
            --bfile $params.genotypes \\
            --keep pqtl_individuals.txt \\
            --freq \\
            --allow-extra-chr

        Rscript $workflow.projectDir/pqtl_conditional_effects_finemap_prepare_loci.R
        """
}

process RUN_COJO {

    errorStrategy "finish"

    label "cojo"

    input:
        path(pqtl_individuals)
        tuple(val(locus), path("locus_files/*"))
    
    output:
        path("${locus}_conditional_signals.txt"), emit: conditional_signal_stats

    script:

        """
        cat locus_files/*.cond.snps | sort | uniq > signal_snps.txt

        plink \\
            --bfile $params.genotypes \\
            --keep $pqtl_individuals \\
            --extract locus_files/${locus}.snps \\
            --r2 \\
            --ld-snp-list signal_snps.txt \\
            --ld-window-kb 2000 \\
            --ld-window 99999 \\
            --ld-window-r2 0.9 \\
            --allow-extra-chr

        awk 'NR > 1 { print \$6; }' plink.ld | sort | uniq > signal_ld_snps.txt

        grep -vwFf signal_ld_snps.txt locus_files/${locus}.snps > non_signal_snps.txt

        plink \\
            --bfile $params.genotypes \\
            --keep $pqtl_individuals \\
            --extract non_signal_snps.txt \\
            --indep 99999 1 100 \\
            --allow-extra-chr
        
        mv signal_snps.txt locus_files/${locus}.snps.pruned
        grep -wFf plink.prune.in non_signal_snps.txt >> locus_files/${locus}.snps.pruned

        for file_name in locus_files/*.cond.snps
        do

            gcta64 \\
                --bfile $params.genotypes \\
                --keep $pqtl_individuals \\
                --extract locus_files/${locus}.snps.pruned \\
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

        cat *.part > ${locus}_conditional_signals.txt
        """
}

process MERGE_CONDITIONAL_OUTPUT {

    label "simple_bash"

    publishDir "$params.output_dir/", mode: "move"

    input:
        path(cis_single_signal_stats)
        path(trans_single_signal_stats)
        path("cis_conditional_signals/*")
        path("trans_conditional_signals/*")
    
    output:
        path("cis_conditional_FINEMAP_CS_summary_statistics.tsv")
        path("trans_conditional_FINEMAP_CS_summary_statistics.tsv")
    
    script:

        """
        # Cis-pQTL
        echo -e "Gene\tSignal\tChr\tSNP\tPosition\tBeta\tSE\tP_Value" > cis_conditional_FINEMAP_CS_summary_statistics.tsv

        cat $cis_single_signal_stats cis_conditional_signals/* >> cis_conditional_FINEMAP_CS_summary_statistics.tsv

        # Trans-pQTL
        echo -e "Gene\tSignal\tChr\tSNP\tPosition\tBeta\tSE\tP_Value" > trans_conditional_FINEMAP_CS_summary_statistics.tsv

        cat $trans_single_signal_stats trans_conditional_signals/* >> trans_conditional_FINEMAP_CS_summary_statistics.tsv
        """
}


workflow {

    PREPARE_LOCI()

    locus_files = PREPARE_LOCI.out.cis_locus_files
        .concat(PREPARE_LOCI.out.trans_locus_files)
        .flatten()
        .map{ file -> [file.getSimpleName().replaceAll(/-\d+/, ""), file] }
        .groupTuple()
    
    RUN_COJO(PREPARE_LOCI.out.mapping_individuals, locus_files)

    conditional_signal_stats = RUN_COJO.out.conditional_signal_stats
        .branch{
            cis: it.getSimpleName() ==~ /ENSG.*/
            trans: it.getSimpleName() ==~ /Trans.*/
        }
    
    MERGE_CONDITIONAL_OUTPUT(
        PREPARE_LOCI.out.cis_single_signal_stats,
        PREPARE_LOCI.out.trans_single_signal_stats,
        conditional_signal_stats.cis.collect(),
        conditional_signal_stats.trans.collect()
    )
}
