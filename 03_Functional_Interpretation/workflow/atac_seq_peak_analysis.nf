nextflow.enable.dsl = 2

params.atlas = "combined"
params.atac_seq_dir = "/nfs/users/nfs_n/nm18/gains_team282/epigenetics/accessibility/merged/atac_seq/"
params.metadata = "/nfs/users/nfs_n/nm18/eQTL_pQTL_Characterization/03_Functional_Interpretation/metadata/reads_atac_seq.txt"
params.peak_saf = "/nfs/users/nfs_n/nm18/gains_team282/epigenetics/accessibility/analysis/atac_seq/combined/consensus_peaks.saf"
params.output_dir = "/nfs/users/nfs_n/nm18/gains_team282/epigenetics/accessibility/analysis/atac_seq/combined/"


//-----------------------------------------------------
// Count Fragments in Peaks
//-----------------------------------------------------

process COUNT_FRAGMENTS {

    errorStrategy "retry"
    maxRetries 5

    label "featureCounts"

    input:
        val(sample)
    
    output:
        path("${sample}.tsv"), emit: peak_fragment_count
        path("${sample}.tsv.summary"), emit: peak_fragment_count_summary
        path("${sample}_macs2.tsv.summary"), emit: macs2_peak_fragment_count_summary
    
    script:

        """
        featureCounts \\
            -p \\
            -a $params.peak_saf -F SAF \\
            -T $task.cpus \\
            -o ${sample}.tsv \\
            $params.atac_seq_dir/$sample/alignment/${sample}.sortedByName.bam
        
        echo -e "GeneID\tChr\tStart\tEnd\tStrand" > macs2_peaks.saf
        awk 'OFS="\t" { print \$1 ":" \$2 "-" \$3, \$1, \$2, \$3, "+"; }' $params.atac_seq_dir/$sample/peaks/${sample}_peaks.narrowPeak >> macs2_peaks.saf

        featureCounts \\
            -p \\
            -a macs2_peaks.saf -F SAF \\
            -T $task.cpus \\
            -o ${sample}_macs2.tsv \\
            $params.atac_seq_dir/$sample/alignment/${sample}.sortedByName.bam
        """
}

process AGGREGATE_FRAGMENT_COUNTS {

    errorStrategy "retry"
    maxRetries 3

    label "simple_bash"

    publishDir "$params.output_dir/", mode: "copy"

    input:
        path("counts/*.tsv")
        path("summaries/*.summary")
        path("macs2_summaries/*.summary")

    output:
        path("peak_counts.tsv")
        path("peak_frips.tsv")
        path("macs2_peak_frips.tsv")

    script:

        """
        for counts_file in counts/*.tsv
        do
            awk 'NR>1 { print \$7; }' \$counts_file > \${counts_file}.part
            sed -i "1s:.*/::" \${counts_file}.part
        done

        paste $params.peak_saf counts/*.part > peak_counts.tsv

        chmod 444 peak_counts.tsv

        echo -e "Sample_File\tAssigned_Reads\tTotal_Reads\tFRiP" > peak_frips.tsv

        for summaries_file in summaries/*.summary
        do
            sample_file=\$(awk 'NR==1 { print \$2; }' \$summaries_file | xargs basename)
            total_reads=\$(awk 'NR>1 { print \$2; }' \$summaries_file | paste -sd+ | bc)
            assigned_reads=\$(awk 'NR==2 { print \$2; }' \$summaries_file)
            frip=\$(echo "\$assigned_reads / \$total_reads" | bc -l)
            echo -e "\$sample_file\t\$assigned_reads\t\$total_reads\t\$frip" > \${summaries_file}.part
        done

        cat summaries/*.part >> peak_frips.tsv

        chmod 444 peak_frips.tsv

        echo -e "Sample_File\tAssigned_Reads\tTotal_Reads\tFRiP" > macs2_peak_frips.tsv

        for summaries_file in macs2_summaries/*.summary
        do
            sample_file=\$(awk 'NR==1 { print \$2; }' \$summaries_file | xargs basename)
            total_reads=\$(awk 'NR>1 { print \$2; }' \$summaries_file | paste -sd+ | bc)
            assigned_reads=\$(awk 'NR==2 { print \$2; }' \$summaries_file)
            frip=\$(echo "\$assigned_reads / \$total_reads" | bc -l)
            echo -e "\$sample_file\t\$assigned_reads\t\$total_reads\t\$frip" > \${summaries_file}.part
        done

        cat macs2_summaries/*.part >> macs2_peak_frips.tsv

        chmod 444 macs2_peak_frips.tsv
        """
}

//-----------------------------------------------------
// Workflow
//-----------------------------------------------------

workflow {

    //-----------------------------------------------------
    // Count Fragments in Peaks
    //-----------------------------------------------------

    samples_list = Channel
        .fromPath(params.metadata)
        .splitCsv(skip: 1)
        .map{row -> row[1]}

    if (params.atlas == "immune") {
        samples_list = samples_list
            .filter{ it =~ /Corces/ || it =~ /Calderon/ }
    }

    if (params.atlas == "neutrophil") {
        samples_list = samples_list
            .filter{ it =~ /Ram-Mohan/ }
    }

    COUNT_FRAGMENTS(samples_list)

    AGGREGATE_FRAGMENT_COUNTS(
        COUNT_FRAGMENTS.out.peak_fragment_count.collect(),
        COUNT_FRAGMENTS.out.peak_fragment_count_summary.collect(),
        COUNT_FRAGMENTS.out.macs2_peak_fragment_count_summary.collect()
    )
}
