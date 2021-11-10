nextflow.enable.dsl = 2

params.atac_seq_dir = "/nfs/users/nfs_n/nm18/gains_team282/epigenetics/calderon_et_al/processed/atac_seq/"
params.output_dir = "/nfs/users/nfs_n/nm18/gains_team282/epigenetics/calderon_et_al/analysis/atac_seq/"

process IDENTIFY_SAMPLES {

    errorStrategy "retry"
    maxRetries 3

    label "simple_bash"

    output:
        path("samples.txt"), emit: samples_file_jaccard
        path("samples.txt"), emit: samples_file_count
    
    script:

    """
    ls -1a $params.atac_seq_dir | grep "SRR" > samples.txt
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
        """
}

//-----------------------------------------------------
// Generate Consensus Peak Set
//-----------------------------------------------------

process CONSENSUS_PEAK_SET {

    errorStrategy "retry"
    maxRetries 3

    label "bedtools"

    publishDir "$params.output_dir/", mode: "copy"

    input:
        path("narrow_peaks/*.narrowPeak")

    output:
        path("consensus_peaks.bed"), emit: consensus_peaks
    
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
        """
}

//-----------------------------------------------------
// Count Fragments in Peaks
//-----------------------------------------------------

process COUNT_FRAGMENTS {

    errorStrategy "retry"
    maxRetries 3

    label "featureCounts"

    publishDir "$params.output_dir/", mode: "copy"

    input:
        path(samples_file)
        path(consensus_peaks)
    
    output:
        path("peak_counts.tsv")
    
    script:

        """
        awk '{ print "$params.atac_seq_dir/" \$1 "/alignment_post_qc/" \$1 ".filtered.bam"; }' $samples_file > alignment_files.txt

        echo -e "GeneID\tChr\tStart\tEnd\tStrand" > consensus_peaks.saf
        awk 'OFS="\t" { print \$1 ":" \$2 "-" \$3, \$1, \$2, \$3, "+"; }' $consensus_peaks >> consensus_peaks.saf

        cat alignment_files.txt | xargs featureCounts -p -a consensus_peaks.saf -F SAF -T $task.cpus -o peak_counts.tsv
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

    SAMPLE_PAIRWISE_JACCARD(IDENTIFY_SAMPLES.out.samples_file_jaccard)

    CONSENSUS_PEAK_SET(peak_files)

    COUNT_FRAGMENTS(IDENTIFY_SAMPLES.out.samples_file_count, CONSENSUS_PEAK_SET.out.consensus_peaks)
}
