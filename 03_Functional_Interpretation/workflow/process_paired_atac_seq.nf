nextflow.enable.dsl = 2

params.output = "/nfs/users/nfs_n/nm18/gains_team282/epigenetics/accessibility/processed/atac_seq/"
params.bowtie2_genome_dir = "/nfs/users/nfs_n/nm18/gains_team282/epigenetics/bowtie2_genome_index/"
params.genome_fasta = "/lustre/scratch118/humgen/resources/rna_seq_genomes/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
params.genome_black_list = "/nfs/users/nfs_n/nm18/eQTL_pQTL_Characterization/03_Functional_Interpretation/data/hg38-blacklist.v2.bed"


//-----------------------------------------------------
// Alignment (A)
//-----------------------------------------------------

process A_BOWTIE2 {

    errorStrategy "retry"
    maxRetries 3

    label "bowtie2"

    input:
        val(name)
        path(fastq_file_1)
        path(fastq_file_2)

    output:
        val(name),                      emit: name
        path("${name}.aligned.bam"),    emit: bam_file

    script:

        def cores = (task.cpus as int) - 1

        """
        bowtie2 \\
            -x $params.bowtie2_genome_dir/Homo_sapiens.GRCh38.99 \\
            -1 $fastq_file_1 -2 $fastq_file_2 \\
            --no-mixed \\
            --no-discordant \\
            --threads $cores \\
            --met-file ${name}.bowtie2.log \\
            | samtools view --bam > ${name}.aligned.bam
        """
}

//-----------------------------------------------------
// Post-Alignment QC (POSTQC)
//-----------------------------------------------------

process POSTQC_SORT {

    errorStrategy "retry"
    maxRetries 3

    label "samtools"

    input:
        val(name)
        path(bam_file)

    output:
        val(name),                                          emit: name
        path("${name}.bam"),                                emit: bam_file
        path("${name}.bam.bai"),                            emit: bam_index_file
        path("*.{bowtie2.log,stats,idxstats,flagstat}"),    emit: mqc_alignment

    script:

        """
        samtools sort -m 1G -@ $task.cpus $bam_file > ${name}.bam
        samtools index -@ $task.cpus ${name}.bam
        samtools idxstats -@ $task.cpus ${name}.bam > ${name}.bam.idxstats
        samtools flagstat -@ $task.cpus ${name}.bam > ${name}.bam.flagstat
        samtools stats -@ $task.cpus ${name}.bam > ${name}.bam.stats
        """
}

process POSTQC_PICARD_METRICS {

    errorStrategy "retry"
    maxRetries 3

    label "picard"

    input:
        val(name)
        path(bam_file)

    output:
        path("*_metrics"), emit: mqc_picard_metrics

    script:

        """
        picard -Xmx${task.memory.toGiga()}g CollectMultipleMetrics \\
            INPUT=$bam_file \\
            OUTPUT=${name}.CollectMultipleMetrics \\
            REFERENCE_SEQUENCE=$params.genome_fasta \\
            VALIDATION_STRINGENCY=LENIENT \\
            TMP_DIR=./tmp/
        """
}

process POSTQC_PICARD_MARK_DUPLICATES {

    errorStrategy "retry"
    maxRetries 3

    label "picard"

    input:
        val(name)
        path(bam_file)

    output:
        val(name),                          emit: name
        path("${name}.duplicates.bam"),     emit: bam_file
        path("${name}.duplicates.bam.bai"), emit: bam_index_file
        path("${name}.MarkDuplicates.txt"), emit: mqc_picard_mark_duplicates

    script:

        """
        picard -Xmx${task.memory.toGiga() - 1}g MarkDuplicates \\
            INPUT=$bam_file \\
            OUTPUT=${name}.duplicates.bam \\
            ASSUME_SORTED=true \\
            REMOVE_DUPLICATES=false \\
            METRICS_FILE=${name}.MarkDuplicates.txt \\
            VALIDATION_STRINGENCY=LENIENT \\
            TMP_DIR=./tmp/

        samtools index -@ $task.cpus ${name}.duplicates.bam
        """
}

process POSTQC_FILTERING {

    errorStrategy "retry"
    maxRetries 3

    label "samtools"

    publishDir "$params.output/$name/alignment_post_qc/", mode: "copy", pattern: "*.{bam,bam.bai}"

    input:
        val(name)
        path(bam_file)
        path(bam_index_file)

    output:
        val(name),                                  emit: name
        path("${name}.filtered.bam"),               emit: bam_file
        path("${name}.filtered.bam.bai"),           emit: bam_index_file
        path("${name}.filtered.sortedByName.bam")
        path("*.{stats,idxstats,flagstat}"),        emit: mqc_post_filtering

    script:

        // --require-flags 0x0001 - Keep reads that are paired.
        // --exclude-flags 0x0004 - Filter out reads that are unmapped.
        // --exclude-flags 0x0200 -F 0x0400 - Filter out reads that are marked as duplicates.
        // --remove-tag XS - XS is the score of the second best alignment, which implies that the read mapped to multiple locations

        """
        samtools view $bam_file \\
            --require-flags 0x0001 \\
            --exclude-flags 0x0004 \\
            --exclude-flags 0x0200 --exclude-flags 0x0400 \\
            --remove-tag XS \\
            --bam \\
            -@ $task.cpus > ${name}.filtered.bam
        
        samtools index -@ $task.cpus ${name}.filtered.bam

        # Filter reads that map to mitochondrial genome (MT)
        samtools idxstats ${name}.filtered.bam | cut -f 1 | grep -v MT | xargs samtools view -b ${name}.filtered.bam > ${name}.filtered.noMT.bam

        # Filter reads in blacklisted region from ENCODE
        bedtools intersect -abam ${name}.filtered.noMT.bam -b $params.genome_black_list -v > ${name}.filtered.bam

        samtools index -@ $task.cpus ${name}.filtered.bam
        samtools idxstats -@ $task.cpus ${name}.filtered.bam > ${name}.filtered.bam.idxstats
        samtools flagstat -@ $task.cpus ${name}.filtered.bam > ${name}.filtered.bam.flagstat
        samtools stats -@ $task.cpus ${name}.filtered.bam > ${name}.filtered.bam.stats

        # Sort reads by name, making it easier to count using featureCounts
        samtools sort -m 1G -@ $task.cpus -n ${name}.filtered.bam > ${name}.filtered.sortedByName.bam
        """
}

//-----------------------------------------------------
// Peak Calling (PEAK)
//-----------------------------------------------------

process PEAK_MACS2_NARROW {

    errorStrategy "retry"
    maxRetries 3

    label "macs2"

    publishDir "$params.output/$name/peaks/", mode: "move"

    input:
        val(name)
        path(bam_file)
        path(bam_index_file)

    output:
        path("${name}_peaks.narrowPeak")
        path("${name}_peaks.tsv")
        path("${name}_summits.bed")

    script:

        """
        macs2 callpeak \\
            -t $bam_file \\
            -f BAMPE \\
            -n $name \\
            --keep-dup all \\
            --nomodel \\
            --nolambda

        mv ${name}_peaks.xls ${name}_peaks.tsv
        """
}

//-----------------------------------------------------
// Summarization (SUMMARY)
//-----------------------------------------------------

process SUMMARY_MULTI_QC {

    errorStrategy "retry"
    maxRetries 3

    label "multiqc"

    publishDir "$params.output/multiqc/alignment/", mode: "move"

    input:
        path("alignment/*")
        path("picard/metrics/*")
        path("picard/metrics/*")
        path("filtering/*")

    output:
        path("*multiqc_report.html")
        path("*_data")

    script:

        """
        cp $workflow.projectDir/multiqc_config.yaml ./

        multiqc -f -d .
        """
}


//-----------------------------------------------------
// Main Workflow
//-----------------------------------------------------

workflow {

    //-----------------------------------------------------
    // Setup
    //-----------------------------------------------------

    align = Channel
        .fromFilePairs("$params.output/*/pre_processing/*_{1,2}.trimmed.fastq.gz")
        .multiMap { name, files ->
            name: name
            fastq_1: files[0]
            fastq_2: files[1]
        }

    //-----------------------------------------------------
    // Alignment (A)
    //-----------------------------------------------------
    
    A_BOWTIE2(align.name, align.fastq_1, align.fastq_2)

    //-----------------------------------------------------
    // Post-Alignment QC (POSTQC)
    //-----------------------------------------------------

    // Sort alignment
    POSTQC_SORT(A_BOWTIE2.out.name, A_BOWTIE2.out.bam_file)

    // Split output from alignment into two channels
    post_align_names = POSTQC_SORT.out.name
        .multiMap { name -> metrics: duplicates: name }
    
    post_align_bam_files = POSTQC_SORT.out.bam_file
        .multiMap { bam_file -> metrics: duplicates: bam_file }

    // Post-Alignment Metrics
    POSTQC_PICARD_METRICS(post_align_names.metrics, post_align_bam_files.metrics)

    // Mark Duplicates
    POSTQC_PICARD_MARK_DUPLICATES(post_align_names.duplicates, post_align_bam_files.duplicates)

    // Filter Reads
    POSTQC_FILTERING(POSTQC_PICARD_MARK_DUPLICATES.out.name, POSTQC_PICARD_MARK_DUPLICATES.out.bam_file, POSTQC_PICARD_MARK_DUPLICATES.out.bam_index_file)

    //-----------------------------------------------------
    // Peak Calling (PEAK)
    //-----------------------------------------------------

    // Call Peaks
    PEAK_MACS2_NARROW(POSTQC_FILTERING.out.name, POSTQC_FILTERING.out.bam_file, POSTQC_FILTERING.out.bam_index_file)

    //-----------------------------------------------------
    // Summarization (SUMMARY)
    //-----------------------------------------------------

    SUMMARY_MULTI_QC(
        POSTQC_SORT.out.mqc_alignment.collect(),
        POSTQC_PICARD_METRICS.out.mqc_picard_metrics.collect(),
        POSTQC_PICARD_MARK_DUPLICATES.out.mqc_picard_mark_duplicates.collect(),
        POSTQC_FILTERING.out.mqc_post_filtering.collect()
    )
}
