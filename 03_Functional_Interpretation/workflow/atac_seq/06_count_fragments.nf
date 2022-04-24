nextflow.enable.dsl = 2

params.atlas = "combined"
params.atac_seq_dir = "/nfs/users/nfs_n/nm18/gains_team282/epigenetics/accessibility/merged/atac_seq/"
params.metadata = "/nfs/users/nfs_n/nm18/eQTL_pQTL_Characterization/03_Functional_Interpretation/metadata/reads_atac_seq.txt"
params.peak_saf = "/nfs/users/nfs_n/nm18/gains_team282/epigenetics/accessibility/analysis/atac_seq/combined/consensus_peaks.saf"
params.peak_set_dir = "/nfs/users/nfs_n/nm18/gains_team282/epigenetics/accessibility/analysis/atac_seq/combined/peak_sets/"
params.cell_type_set_dir = "/nfs/users/nfs_n/nm18/gains_team282/epigenetics/accessibility/analysis/atac_seq/cell_type_peak_sets/"
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
        val(group)
        val(cell_type)
    
    output:
        path("${sample}.tsv"),                                              emit: peak_fragment_count
        path("${sample}.tsv.summary"),                                      emit: peak_fragment_count_summary
        tuple val(group), path("${sample}.${group}.tsv"),                   emit: group_peak_fragment_count
        tuple val(group), path("${sample}.${group}.tsv.summary"),           emit: group_peak_fragment_count_summary
        tuple val(cell_type), path("${sample}.${cell_type}.tsv"),           emit: cell_type_peak_fragment_count
        tuple val(cell_type), path("${sample}.${cell_type}.tsv.summary"),   emit: cell_type_peak_fragment_count_summary
        path("${sample}_macs2.tsv.summary"),                                emit: macs2_peak_fragment_count_summary
    
    script:

        // featureCounts arguments
        // -p - Paired End Reads
        // -O - Count Reads that overlap Multiple Regions

        """
        featureCounts \\
            -p \\
            -O \\
            -a $params.peak_saf -F SAF \\
            -T $task.cpus \\
            -o ${sample}.tsv \\
            $params.atac_seq_dir/$sample/alignment/${sample}.sortedByName.bam

        if [ -s $params.peak_set_dir/${group}.peaks.bed ]
        then

            featureCounts \\
                -p \\
                -O \\
                -a $params.peak_set_dir/${group}.peaks.saf -F SAF \\
                -T $task.cpus \\
                -o ${sample}.${group}.tsv \\
                $params.atac_seq_dir/$sample/alignment/${sample}.sortedByName.bam

        else

            # It is possible that no peaks passed the q-value cutoff for this biological group
            echo "# Program:featureCounts" > ${sample}.${group}.tsv
            echo -e "Geneid\\tChr\\tStart\\tEnd\\tStrand\\tLength\\t$params.atac_seq_dir/$sample/alignment/${sample}.sortedByName.bam" >> ${sample}.${group}.tsv

            num_reads=\$(samtools view -c $params.atac_seq_dir/$sample/alignment/${sample}.sortedByName.bam)

            echo -e "Status\\t$params.atac_seq_dir/$sample/alignment/${sample}.sortedByName.bam" > ${sample}.${group}.tsv.summary
            echo -e "Assigned\\t0" >> ${sample}.${group}.tsv.summary
            echo -e "Unassigned_Unmapped\\t0" >> ${sample}.${group}.tsv.summary
            echo -e "Unassigned_Read_Type\\t0" >> ${sample}.${group}.tsv.summary
            echo -e "Unassigned_Singleton\\t0" >> ${sample}.${group}.tsv.summary
            echo -e "Unassigned_MappingQuality\\t0" >> ${sample}.${group}.tsv.summary
            echo -e "Unassigned_Chimera\\t0" >> ${sample}.${group}.tsv.summary
            echo -e "Unassigned_FragmentLength\\t0" >> ${sample}.${group}.tsv.summary
            echo -e "Unassigned_Duplicate\\t0" >> ${sample}.${group}.tsv.summary
            echo -e "Unassigned_MultiMapping\\t0" >> ${sample}.${group}.tsv.summary
            echo -e "Unassigned_Secondary\\t0" >> ${sample}.${group}.tsv.summary
            echo -e "Unassigned_NonSplit\\t0" >> ${sample}.${group}.tsv.summary
            echo -e "Unassigned_NoFeatures\\t\${num_reads}" >> ${sample}.${group}.tsv.summary
            echo -e "Unassigned_Overlapping_Length\\t0" >> ${sample}.${group}.tsv.summary
            echo -e "Unassigned_Ambiguity\\t0" >> ${sample}.${group}.tsv.summary

        fi

        featureCounts \\
            -p \\
            -O \\
            -a $params.cell_type_set_dir/${cell_type}.peaks.saf -F SAF \\
            -T $task.cpus \\
            -o ${sample}.${cell_type}.tsv \\
            $params.atac_seq_dir/$sample/alignment/${sample}.sortedByName.bam
        
        echo -e "GeneID\tChr\tStart\tEnd\tStrand" > macs2_peaks.saf
        awk 'OFS="\t" { print \$1 ":" \$2 "-" \$3, \$1, \$2, \$3, "+"; }' $params.atac_seq_dir/$sample/peaks/${sample}_peaks.narrowPeak >> macs2_peaks.saf

        featureCounts \\
            -p \\
            -O \\
            -a macs2_peaks.saf -F SAF \\
            -T $task.cpus \\
            -o ${sample}_macs2.tsv \\
            $params.atac_seq_dir/$sample/alignment/${sample}.sortedByName.bam
        """
}

process AGGREGATE_GROUP_FRAGMENT_COUNTS {

    errorStrategy "retry"
    maxRetries 3

    label "simple_bash"

    publishDir "$params.output_dir/peak_sets/", mode: "copy"

    input:
        val(group)
        path("counts/*")
        path("summaries/*")

    output:
        path("${group}_peak_counts.tsv")
        path("${group}_peak_frips.tsv")

    script:

        """
        for counts_file in counts/*.tsv
        do
            awk 'NR>1 { print \$7; }' \$counts_file > \${counts_file}.part
            sed -i "1s:.*/::" \${counts_file}.part
            sed -i 's/.sortedByName.bam//g' \${counts_file}.part
        done

        paste $params.peak_set_dir/${group}.peaks.saf counts/*.part > ${group}_peak_counts.tsv

        echo -e "Sample_File\tAssigned_Reads\tTotal_Reads\tFRiP" > ${group}_peak_frips.tsv

        for summaries_file in summaries/*.summary
        do
            sample_file=\$(awk 'NR==1 { print \$2; }' \$summaries_file | xargs basename)
            total_reads=\$(awk 'NR>1 { print \$2; }' \$summaries_file | paste -sd+ | bc)
            assigned_reads=\$(awk 'NR==2 { print \$2; }' \$summaries_file)
            frip=\$(echo "\$assigned_reads / \$total_reads" | bc -l)
            echo -e "\$sample_file\t\$assigned_reads\t\$total_reads\t\$frip" > \${summaries_file}.part
            sed -i 's/.sortedByName.bam//g' \${summaries_file}.part
        done

        cat summaries/*.part >> ${group}_peak_frips.tsv
        """
}

process AGGREGATE_CELL_TYPE_FRAGMENT_COUNTS {

    errorStrategy "retry"
    maxRetries 3

    label "simple_bash"

    publishDir "$params.output_dir/cell_type_peak_sets/", mode: "copy"

    input:
        val(cell_type)
        path("counts/*")
        path("summaries/*")

    output:
        path("${cell_type}_peak_counts.tsv")
        path("${cell_type}_peak_frips.tsv")

    script:

        """
        for counts_file in counts/*.tsv
        do
            awk 'NR>1 { print \$7; }' \$counts_file > \${counts_file}.part
            sed -i "1s:.*/::" \${counts_file}.part
            sed -i 's/.sortedByName.bam//g' \${counts_file}.part
        done

        paste $params.cell_type_set_dir/${cell_type}.peaks.saf counts/*.part > ${cell_type}_peak_counts.tsv

        echo -e "Sample_File\tAssigned_Reads\tTotal_Reads\tFRiP" > ${cell_type}_peak_frips.tsv

        for summaries_file in summaries/*.summary
        do
            sample_file=\$(awk 'NR==1 { print \$2; }' \$summaries_file | xargs basename)
            total_reads=\$(awk 'NR>1 { print \$2; }' \$summaries_file | paste -sd+ | bc)
            assigned_reads=\$(awk 'NR==2 { print \$2; }' \$summaries_file)
            frip=\$(echo "\$assigned_reads / \$total_reads" | bc -l)
            echo -e "\$sample_file\t\$assigned_reads\t\$total_reads\t\$frip" > \${summaries_file}.part
            sed -i 's/.sortedByName.bam//g' \${summaries_file}.part
        done

        cat summaries/*.part >> ${cell_type}_peak_frips.tsv
        """
}

process AGGREGATE_CONSENSUS_FRAGMENT_COUNTS {

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
            sed -i 's/.sortedByName.bam//g' \${counts_file}.part
        done

        paste $params.peak_saf counts/*.part > peak_counts.tsv

        echo -e "Sample_File\tAssigned_Reads\tTotal_Reads\tFRiP" > peak_frips.tsv

        for summaries_file in summaries/*.summary
        do
            sample_file=\$(awk 'NR==1 { print \$2; }' \$summaries_file | xargs basename)
            total_reads=\$(awk 'NR>1 { print \$2; }' \$summaries_file | paste -sd+ | bc)
            assigned_reads=\$(awk 'NR==2 { print \$2; }' \$summaries_file)
            frip=\$(echo "\$assigned_reads / \$total_reads" | bc -l)
            echo -e "\$sample_file\t\$assigned_reads\t\$total_reads\t\$frip" > \${summaries_file}.part
            sed -i 's/.sortedByName.bam//g' \${summaries_file}.part
        done

        cat summaries/*.part >> peak_frips.tsv

        echo -e "Sample_File\tAssigned_Reads\tTotal_Reads\tFRiP" > macs2_peak_frips.tsv

        for summaries_file in macs2_summaries/*.summary
        do
            sample_file=\$(awk 'NR==1 { print \$2; }' \$summaries_file | xargs basename)
            total_reads=\$(awk 'NR>1 { print \$2; }' \$summaries_file | paste -sd+ | bc)
            assigned_reads=\$(awk 'NR==2 { print \$2; }' \$summaries_file)
            frip=\$(echo "\$assigned_reads / \$total_reads" | bc -l)
            echo -e "\$sample_file\t\$assigned_reads\t\$total_reads\t\$frip" > \${summaries_file}.part
            sed -i 's/.sortedByName.bam//g' \${summaries_file}.part
        done

        cat macs2_summaries/*.part >> macs2_peak_frips.tsv
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
        .unique{ row -> row[1] }

    if (params.atlas == "immune") {
        samples_list = samples_list
            .filter{ row -> row[3] =~ /Corces/ || row[3] =~ /Calderon/ }
            .multiMap{ row ->
                sample: row[1]
                group: row[5] + (row[6] == "TRUE" ? "-Stimulated" : "-Control")
                cell_type: row[5].replaceAll(/-.*/, "")
            }
    }

    if (params.atlas == "neutrophil") {
        samples_list = samples_list
            .filter{ row -> row[3] =~ /Ram-Mohan/ }
            .multiMap{ row ->
                sample: row[1]
                group: row[5] + "-" + row[7]
                cell_type: row[5].replaceAll(/-.*/, "")
            }
    }

    COUNT_FRAGMENTS(samples_list.sample, samples_list.group, samples_list.cell_type)

    group_counts = COUNT_FRAGMENTS.out.group_peak_fragment_count
        .join(COUNT_FRAGMENTS.out.group_peak_fragment_count_summary)
        .groupTuple()
        .multiMap{ item ->
            group: item[0]
            counts: item[1]
            summaries: item[2]
        }
    
    AGGREGATE_GROUP_FRAGMENT_COUNTS(
        group_counts.group,
        group_counts.counts,
        group_counts.summaries
    )

    cell_type_counts = COUNT_FRAGMENTS.out.cell_type_peak_fragment_count
        .join(COUNT_FRAGMENTS.out.cell_type_peak_fragment_count_summary)
        .groupTuple()
        .multiMap{ item ->
            cell_type: item[0]
            counts: item[1]
            summaries: item[2]
        }

    AGGREGATE_CELL_TYPE_FRAGMENT_COUNTS(
        cell_type_counts.cell_type,
        cell_type_counts.counts,
        cell_type_counts.summaries
    )

    AGGREGATE_CONSENSUS_FRAGMENT_COUNTS(
        COUNT_FRAGMENTS.out.peak_fragment_count.collect(),
        COUNT_FRAGMENTS.out.peak_fragment_count_summary.collect(),
        COUNT_FRAGMENTS.out.macs2_peak_fragment_count_summary.collect()
    )
}
