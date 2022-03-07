nextflow.enable.dsl = 2

params.atac_seq_dir = "/nfs/users/nfs_n/nm18/gains_team282/epigenetics/accessibility/merged/atac_seq/"
params.metadata = "/nfs/users/nfs_n/nm18/eQTL_pQTL_Characterization/03_Functional_Interpretation/metadata/reads_atac_seq.txt"
params.genome_chr_lengths = "/nfs/users/nfs_n/nm18/gains_team282/epigenetics/star_genome_index/chrNameLength.txt"
params.genome_annotation = "/lustre/scratch118/humgen/resources/rna_seq_genomes/Homo_sapiens.GRCh38.99.gtf"
params.output_dir = "/nfs/users/nfs_n/nm18/gains_team282/epigenetics/accessibility/analysis/atac_seq/"


//-----------------------------------------------------
// Sample Peak QC
//-----------------------------------------------------

process SAMPLE_PEAK_WIDTHS {

    errorStrategy "retry"
    maxRetries 3

    label "simple_bash"

    publishDir "$params.output_dir/", mode: "copy"

    output:
        path("peak_widths.tsv")

    script:
        
        // Awk quotations do not play well with GNU parallel.
        // I have wrapped the awk command into a bash function.
        """
        awk -F ',' 'NR > 1 { print \$2; }' $params.metadata | sort | uniq > samples.txt

        awk_command() {
            SAMPLE=\$1
            awk -v sample=\$SAMPLE 'OFS="\t" { print sample, \$3 - \$2; }' $params.atac_seq_dir/\${SAMPLE}/peaks/\${SAMPLE}_peaks.narrowPeak > \${SAMPLE}.part
        }
        export -f awk_command

        parallel \\
            -a samples.txt \\
            "awk_command {1}"

        echo -e "Sample\tPeak_Width" > peak_widths.tsv
        cat *.part >> peak_widths.tsv

        chmod 444 peak_widths.tsv
        """
}

process SAMPLE_PAIRWISE_JACCARD {

    errorStrategy "retry"
    maxRetries 3

    label "bedtools"

    publishDir "$params.output_dir/", mode: "copy"

    output:
        path("jaccard_values.tsv")

    script:

        """
        awk -F ',' 'NR > 1 { print \$2; }' $params.metadata | sort | uniq > samples.txt

        parallel \\
            -a samples.txt \\
            -a samples.txt \\
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
// Sample TSS Enrichment Scores
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

        grep "protein_coding" $params.genome_annotation | awk -F '\t' 'OFS="\t" {
                if (\$3 == "gene") {
                    if (\$7 == "+" && \$4 - 2000 > 0) { print \$1, \$4 - 2000, \$4 + 1999, \$1 ":" \$4, "0", "+"; }
                    if (\$7 == "-" && \$5 - 1999 > 0) { print \$1, \$5 - 1999, \$5 + 2000, \$1 ":" \$5, "0", "-"; }
                }
            }' > tss.2000.bed
        bedtools intersect -a tss.2000.bed -b chromosomes.bed -wa -f 1 | uniq -u > tss.2000.filtered.bed

        bedtools slop -i tss.2000.filtered.bed -g $params.genome_chr_lengths -b -1950 > tss.100.filtered.bed

        bedtools slop -i tss.2000.filtered.bed -g $params.genome_chr_lengths -b -99 | \\
            bedtools flank -i - -g $params.genome_chr_lengths -l 99 -r 0 > tss.100.left.flank.filtered.bed

        bedtools slop -i tss.2000.filtered.bed -g $params.genome_chr_lengths -b -99 | \\
            bedtools flank -i - -g $params.genome_chr_lengths -r 99 -l 0 > tss.100.right.flank.filtered.bed

        echo -e "GeneID\tChr\tStart\tEnd\tStrand" > tss.100.saf
        awk -F '\t' 'OFS="\t" {
            print \$4, \$1, \$2, \$3, \$6;
            }' tss.100.filtered.bed >> tss.100.saf

        echo -e "GeneID\tChr\tStart\tEnd\tStrand" > tss.100.left.flank.saf
        awk -F '\t' 'OFS="\t" {
            print \$4, \$1, \$2, \$3, \$6;
            }' tss.100.left.flank.filtered.bed >> tss.100.left.flank.saf
        
        echo -e "GeneID\tChr\tStart\tEnd\tStrand" > tss.100.right.flank.saf
        awk -F '\t' 'OFS="\t" {
            print \$4, \$1, \$2, \$3, \$6;
            }' tss.100.right.flank.filtered.bed >> tss.100.right.flank.saf
        """
}

process TSS_ENRICHMENT_SCORES {

    errorStrategy "retry"
    maxRetries 5

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
            $params.atac_seq_dir/$sample/alignment/${sample}.sortedByName.bam
        
        featureCounts \\
            -p \\
            -a $tss_100_left -F SAF \\
            -T $task.cpus \\
            -o tss_left_flank_counts.tsv \\
            $params.atac_seq_dir/$sample/alignment/${sample}.sortedByName.bam

        featureCounts \\
            -p \\
            -a $tss_100_right -F SAF \\
            -T $task.cpus \\
            -o tss_right_flank_counts.tsv \\
            $params.atac_seq_dir/$sample/alignment/${sample}.sortedByName.bam
        
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
// Consensus Peak Set
//-----------------------------------------------------

process CONSENSUS_PEAK_SET {

    errorStrategy "retry"
    maxRetries 3

    label "bedtools"

    publishDir "$params.output_dir/$atlas/", mode: "copy", pattern: "*.{bed,saf}"

    input:
        val(atlas)
        path("narrow_peaks/*.narrowPeak")

    output:
        tuple val(atlas), path("consensus_peaks.bed"), emit: consensus_peaks_coverage
        tuple val(atlas), path("consensus_peaks.saf")
    
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

        echo -e "GeneID\tChr\tStart\tEnd\tStrand" > consensus_peaks.saf
        awk 'OFS="\t" { print \$1 ":" \$2 "-" \$3, \$1, \$2, \$3, "+"; }' consensus_peaks.bed >> consensus_peaks.saf
        """
}

process COUNT_SAMPLES {

    errorStrategy "retry"
    maxRetries 3

    label "bedtools"

    publishDir "$params.output_dir/$atlas/", mode: "copy"

    input:
        tuple val(atlas), path(consensus_peaks)

    output:
        path("peak_coverage.tsv")
    
    script:

        """
        awk -F ',' 'NR > 1 { print \$2; }' $params.metadata | sort | uniq > samples.txt

        sort -k1,1 -k2,2n $consensus_peaks > consensus_peaks.sorted.bed

        parallel \\
            -a samples.txt \\
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
// Workflow
//-----------------------------------------------------

workflow {

    //-----------------------------------------------------
    // Sample Peak QC
    //-----------------------------------------------------

    SAMPLE_PEAK_WIDTHS()

    SAMPLE_PAIRWISE_JACCARD()

    //-----------------------------------------------------
    // Sample TSS Enrichment Scores
    //-----------------------------------------------------

    samples_list = Channel
        .fromPath(params.metadata)
        .splitCsv(skip: 1)
        .map{row -> row[1]}
        .unique()

    IDENTIFY_TSS_REGIONS()

    TSS_ENRICHMENT_SCORES(
        samples_list,
        IDENTIFY_TSS_REGIONS.out.tss_100,
        IDENTIFY_TSS_REGIONS.out.tss_100_left,
        IDENTIFY_TSS_REGIONS.out.tss_100_right
    )

    AGGREGATE_TSS_ENRICHMENT_SCORES(TSS_ENRICHMENT_SCORES.out.tss_enrichment.collect())

    //-----------------------------------------------------
    // Consensus Peak Set
    //-----------------------------------------------------

    peak_files = Channel
        .from(
            ["neutrophil", file("$params.atac_seq_dir/Ram-Mohan-*/peaks/*.narrowPeak")],
            ["immune", file("$params.atac_seq_dir/{Corces,Calderon}-*/peaks/*.narrowPeak")],
            ["combined", file("$params.atac_seq_dir/**.narrowPeak")]
        )
        .multiMap{ pair -> 
            atlas: pair[0]
            peaks: pair[1]
        }

    CONSENSUS_PEAK_SET(peak_files.atlas, peak_files.peaks)

    COUNT_SAMPLES(CONSENSUS_PEAK_SET.out.consensus_peaks_coverage)
}
