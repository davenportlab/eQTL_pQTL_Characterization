nextflow.enable.dsl = 2

params.atac_seq_dir = "/nfs/users/nfs_n/nm18/gains_team282/epigenetics/calderon_et_al/processed/atac_seq/"
params.genome_chr_lengths = "/nfs/users/nfs_n/nm18/gains_team282/epigenetics/star_genome_index/chrNameLength.txt"
params.genome_annotation = "/lustre/scratch118/humgen/resources/rna_seq_genomes/Homo_sapiens.GRCh38.99.gtf"
params.output_dir = "/nfs/users/nfs_n/nm18/gains_team282/epigenetics/calderon_et_al/analysis/atac_seq/"

process IDENTIFY_SAMPLES {

    errorStrategy "retry"
    maxRetries 3

    label "simple_bash"

    output:
        path("samples.txt"), emit: samples_file_widths
        path("samples.txt"), emit: samples_file_jaccard
        path("samples.txt"), emit: samples_file_coverage
        path("samples.txt"), emit: samples_file_count
    
    script:

    """
    ls -1a $params.atac_seq_dir | grep "SRR" > samples.txt
    """
}

//-----------------------------------------------------
// Sample Peak Widths
//-----------------------------------------------------

process SAMPLE_PEAK_WIDTHS {

    errorStrategy "retry"
    maxRetries 3

    label "simple_bash"

    publishDir "$params.output_dir/", mode: "copy"

    input:
        path(samples_file)

    output:
        path("peak_widths.tsv")

    script:
        
        // Awk quotations do not play well with GNU parallel.
        // I have wrapped the awk command into a bash function.
        """
        awk_command() {
            SAMPLE=\$1
            awk -v sample=\$SAMPLE 'OFS="\t" { print sample, \$3 - \$2; }' $params.atac_seq_dir/\${SAMPLE}/peaks/\${SAMPLE}_peaks.narrowPeak > \${SAMPLE}.part
        }
        export -f awk_command

        parallel \\
            -a $samples_file \\
            "awk_command {1}"

        echo -e "Sample\tPeak_Width" > peak_widths.tsv
        cat *.part >> peak_widths.tsv

        chmod 444 peak_widths.tsv
        """
}

//-----------------------------------------------------
// Comparing Peak Sets from Samples
//-----------------------------------------------------

process SAMPLE_PAIRWISE_JACCARD {

    errorStrategy "retry"
    maxRetries 3

    label "bedtools"

    publishDir "$params.output_dir/", mode: "copy"

    input:
        path(samples_file)

    output:
        path("jaccard_values.tsv")

    script:

        """
        parallel \\
            -a $samples_file \\
            -a $samples_file \\
            "bedtools jaccard \\
            -a $params.atac_seq_dir/{1}/peaks/{1}_peaks.narrowPeak -b $params.atac_seq_dir/{2}/peaks/{2}_peaks.narrowPeak \\
            | awk 'NR > 1' \\
            | cut -f 3 \\
            | sed 's/\$/\t{1}\t{2}/'" \\
            > jaccard_values.tsv

        chmod 444 jaccard_values.tsv
        """
}

//-----------------------------------------------------
// Generate Consensus Peak Set
//-----------------------------------------------------

process CONSENSUS_PEAK_SET {

    errorStrategy "retry"
    maxRetries 3

    label "bedtools"

    publishDir "$params.output_dir/", mode: "copy", pattern: "*.bed"

    input:
        path("narrow_peaks/*.narrowPeak")

    output:
        path("consensus_peaks.bed"), emit: consensus_peaks_coverage
        path("consensus_peaks.saf"), emit: consensus_peaks_count
        path("consensus_peaks.saf"), emit: consensus_peaks_aggregate
    
    script:

        // -d 0 - Maximum distance allowed is 0 / Overlap must be at least 1bp
        // -c 1,7,8,9,10 - Operate on columns 1 (Chromosome Name), 7 (Fold Enrichment), 8 (-log10(P-Value)), 9 (-log10(Q-Value)), and 10 (Point-Source for Peak)
        // -o count,collapse,collapse,collapse,collapse - Count peaks using column 1 and collapse columns 7,8,9,10 into a comma-separated list

        // Filtering via AWK
        //  1. Any peak that is wider than 3 kb is removed
        //  2. Any peak only present in one sample is removed
        """
        cat narrow_peaks/*.narrowPeak > narrow_peaks.bed

        sort -k1,1 -k2,2n narrow_peaks.bed > narrow_peaks.sorted.bed

        bedtools merge \\
            -d 0 \\
            -c 1,7,8,9,10 \\
            -o count,collapse,collapse,collapse,collapse \\
            -i narrow_peaks.sorted.bed | \\
            awk -F '\t' -v OFS='\t' '{
                if (\$3 - \$2 <= 3000) {
                    if (\$4 >= 2) {
                        print \$0;
                    }
                }
            }' > consensus_peaks.bed
        
        chmod 444 consensus_peaks.bed

        echo -e "PeakID\tChr\tStart\tEnd\tStrand" > consensus_peaks.saf
        awk 'OFS="\t" { print \$1 ":" \$2 "-" \$3, \$1, \$2, \$3, "+"; }' consensus_peaks.bed >> consensus_peaks.saf
        """
}

//-----------------------------------------------------
// Count Samples per Consensus Peak
//-----------------------------------------------------

process COUNT_SAMPLES {

    errorStrategy "retry"
    maxRetries 3

    label "bedtools"

    publishDir "$params.output_dir/", mode: "copy"

    input:
        path(samples_file)
        path(consensus_peaks)

    output:
        path("peak_coverage.tsv")
    
    script:

        """
        sort -k1,1 -k2,2n $consensus_peaks > consensus_peaks.sorted.bed

        parallel \\
            -a $samples_file \\
            "echo '{1}' > {1}.part; \\
            sort -k1,1 -k2,2n $params.atac_seq_dir/{1}/peaks/{1}_peaks.narrowPeak \\
            | bedtools coverage -a consensus_peaks.sorted.bed -b stdin -sorted \\
            | cut -f 9 \\
            >> {1}.part"

        echo -e "Chr\tStart\tEnd" > peak_info.tsv
        awk 'OFS="\t" { print \$1, \$2, \$3; }' consensus_peaks.sorted.bed >> peak_info.tsv

        paste peak_info.tsv *.part > peak_coverage.tsv

        chmod 444 peak_coverage.tsv
        """
}

//-----------------------------------------------------
// TSS Enrichment Scores of Samples
//-----------------------------------------------------

process IDENTIFY_TSS_REGIONS {

    errorStrategy "retry"
    maxRetries 3

    label "simple_bash"

    output:
        path("tss.100.saf"),                emit: tss_100
        path("tss.100.left.flank.saf"),     emit: tss_100_left
        path("tss.100.right.flank.saf"),    emit: tss_100_right
    
    script:

        """
        awk -F '\t' 'OFS="\t" { print \$1, "1", \$2; }' $params.genome_chr_lengths > chromosomes.bed

        awk -F '\t' 'OFS="\t" {
                if (\$3 == "gene") {
                    if (\$7 == "+" && \$4 - 2000 > 0) { print \$1, \$4 - 2000, \$4 + 1999, \$1 ":" \$4, "0", "+"; }
                    if (\$7 == "-" && \$5 - 1999 > 0) { print \$1, \$5 - 1999, \$5 + 2000, \$1 ":" \$5, "0", "-"; }
                }
            }' $params.genome_annotation > tss.2000.bed
        bedtools intersect -a tss.2000.bed -b chromosomes.bed -wa -f 1 | uniq -u > tss.2000.filtered.bed

        bedtools slop -i tss.2000.filtered.bed -g $params.genome_chr_lengths -b -1950 > tss.100.filtered.bed

        bedtools slop -i tss.2000.filtered.bed -g $params.genome_chr_lengths -b -99 | \\
            bedtools flank -i - -g $params.genome_chr_lengths -l 99 -r 0 > tss.100.left.flank.filtered.bed

        bedtools slop -i tss.2000.filtered.bed -g $params.genome_chr_lengths -b -99 | \\
            bedtools flank -i - -g $params.genome_chr_lengths -r 99 -l 0 > tss.100.right.flank.filtered.bed

        echo -e "PeakID\tChr\tStart\tEnd\tStrand" > tss.100.saf
        awk -F '\t' 'OFS="\t" {
            print \$4, \$1, \$2, \$3, \$6;
            }' tss.100.filtered.bed >> tss.100.saf

        echo -e "PeakID\tChr\tStart\tEnd\tStrand" > tss.100.left.flank.saf
        awk -F '\t' 'OFS="\t" {
            print \$4, \$1, \$2, \$3, \$6;
            }' tss.100.left.flank.filtered.bed >> tss.100.left.flank.saf
        
        echo -e "PeakID\tChr\tStart\tEnd\tStrand" > tss.100.right.flank.saf
        awk -F '\t' 'OFS="\t" {
            print \$4, \$1, \$2, \$3, \$6;
            }' tss.100.right.flank.filtered.bed >> tss.100.right.flank.saf
        """
}

process TSS_ENRICHMENT_SCORES {

    errorStrategy "retry"
    maxRetries 3

    label "featureCounts"

    input:
        val(sample)
        path(tss_100)
        path(tss_100_left)
        path(tss_100_right)
    
    output:
        path("${sample}_tss_enrichment.tsv"), emit: tss_enrichment

    script:

        """
        featureCounts \\
            -p \\
            -a $tss_100 -F SAF \\
            -T $task.cpus \\
            -o tss_counts.tsv \\
            $params.atac_seq_dir/$sample/alignment_post_qc/${sample}.filtered.sortedByName.bam
        
        featureCounts \\
            -p \\
            -a $tss_100_left -F SAF \\
            -T $task.cpus \\
            -o tss_left_flank_counts.tsv \\
            $params.atac_seq_dir/$sample/alignment_post_qc/${sample}.filtered.sortedByName.bam

        featureCounts \\
            -p \\
            -a $tss_100_right -F SAF \\
            -T $task.cpus \\
            -o tss_right_flank_counts.tsv \\
            $params.atac_seq_dir/$sample/alignment_post_qc/${sample}.filtered.sortedByName.bam
        
        TSS=\$(awk '{ s += \$7; } END { print s; }' tss_counts.tsv)
        LEFT=\$(awk '{ s += \$7; } END { print s; }' tss_left_flank_counts.tsv)
        RIGHT=\$(awk '{ s += \$7; } END { print s; }' tss_right_flank_counts.tsv)

        TSS_ENRICHMENT=\$(echo "(2 * \$TSS) / (\$LEFT + \$RIGHT)" | bc -l) 
        echo -e "$sample\t\$TSS_ENRICHMENT" > ${sample}_tss_enrichment.tsv
        """
}

process AGGREGATE_TSS_ENRICHMENT_SCORES {

    errorStrategy "retry"
    maxRetries 3

    label "simple_bash"

    publishDir "$params.output_dir/", mode: "copy"

    input:
        path("scores/*.tsv")
    
    output:
        path("tss_enrichment_scores.tsv")

    script:

        """
        echo -e "Sample\tTSS_Enrichment_Score" > tss_enrichment_scores.tsv
        cat scores/*.tsv >> tss_enrichment_scores.tsv
        chmod 444 tss_enrichment_scores.tsv
        """
}

//-----------------------------------------------------
// Count Fragments in Peaks
//-----------------------------------------------------

process COUNT_FRAGMENTS {

    errorStrategy "retry"
    maxRetries 3

    label "featureCounts"

    input:
        val(sample)
        path(consensus_peaks)
    
    output:
        path("${sample}.tsv"),                  emit: peak_fragment_count
        path("${sample}.tsv.summary"),          emit: peak_fragment_count_summary
        path("${sample}_macs2.tsv.summary"),    emit: macs2_peak_fragment_count_summary
    
    script:

        """
        featureCounts \\
            -p \\
            -a $consensus_peaks -F SAF \\
            -T $task.cpus \\
            -o ${sample}.tsv \\
            $params.atac_seq_dir/$sample/alignment_post_qc/${sample}.filtered.sortedByName.bam
        
        echo -e "PeakID\tChr\tStart\tEnd\tStrand" > macs2_peaks.saf
        awk 'OFS="\t" { print \$1 ":" \$2 "-" \$3, \$1, \$2, \$3, "+"; }' $params.atac_seq_dir/$sample/peaks/${sample}_peaks.narrowPeak >> macs2_peaks.saf

        featureCounts \\
            -p \\
            -a macs2_peaks.saf -F SAF \\
            -T $task.cpus \\
            -o ${sample}_macs2.tsv \\
            $params.atac_seq_dir/$sample/alignment_post_qc/${sample}.filtered.sortedByName.bam
        """
}

process AGGREGATE_FRAGMENT_COUNTS {

    errorStrategy "retry"
    maxRetries 3

    label "simple_bash"

    publishDir "$params.output_dir/", mode: "copy"

    input:
        path(consensus_peaks)
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

        paste $consensus_peaks counts/*.part > peak_counts.tsv

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

    IDENTIFY_SAMPLES()

    peak_files = Channel
        .fromPath("$params.atac_seq_dir/**.narrowPeak")
        .collect()

    SAMPLE_PEAK_WIDTHS(IDENTIFY_SAMPLES.out.samples_file_widths)

    SAMPLE_PAIRWISE_JACCARD(IDENTIFY_SAMPLES.out.samples_file_jaccard)

    CONSENSUS_PEAK_SET(peak_files)

    COUNT_SAMPLES(IDENTIFY_SAMPLES.out.samples_file_coverage, CONSENSUS_PEAK_SET.out.consensus_peaks_coverage)

    samples_list = IDENTIFY_SAMPLES.out.samples_file_count
        .splitText()
        .map{ sample -> sample.trim() }
        .multiMap{ sample -> tss_enrichment_scores: count_fragment: sample }

    IDENTIFY_TSS_REGIONS()

    TSS_ENRICHMENT_SCORES(
        samples_list.tss_enrichment_scores,
        IDENTIFY_TSS_REGIONS.out.tss_100,
        IDENTIFY_TSS_REGIONS.out.tss_100_left,
        IDENTIFY_TSS_REGIONS.out.tss_100_right
    )

    AGGREGATE_TSS_ENRICHMENT_SCORES(TSS_ENRICHMENT_SCORES.out.tss_enrichment.collect())

    COUNT_FRAGMENTS(samples_list.count_fragment, CONSENSUS_PEAK_SET.out.consensus_peaks_count)

    AGGREGATE_FRAGMENT_COUNTS(
        CONSENSUS_PEAK_SET.out.consensus_peaks_aggregate, 
        COUNT_FRAGMENTS.out.peak_fragment_count.collect(),
        COUNT_FRAGMENTS.out.peak_fragment_count_summary.collect(),
        COUNT_FRAGMENTS.out.macs2_peak_fragment_count_summary.collect()
    )
}
