nextflow.enable.dsl = 2

params.output = "/nfs/users/nfs_n/nm18/gains_team282/epigenetics/calderon_et_al/processed/rna_seq/"
params.star_genome_dir = "/nfs/users/nfs_n/nm18/gains_team282/epigenetics/star_genome_index/"
params.genome_fasta = "/lustre/scratch118/humgen/resources/rna_seq_genomes/Homo_sapiens.GRCh38.dna.primary_assembly.fa"


//-----------------------------------------------------
// Alignment (A)
//-----------------------------------------------------

process A_STAR {

    errorStrategy "retry"
    maxRetries 3

    label "star"

    publishDir "$params.output/$name/alignment/", mode: "copy", pattern: "*.ReadsPerGene.out.tab"

    input:
        val(name)
        path(fastq_file_1)
        path(fastq_file_2)

    output:
        val(name),                                      emit: name_picard_metrics
        val(name),                                      emit: name_picard_mark_duplicates
        path("${name}.Aligned.sortedByCoord.out.bam"),  emit: bam_file_picard_metrics
        path("${name}.Aligned.sortedByCoord.out.bam"),  emit: bam_file_picard_mark_duplicates
        path("*.ReadsPerGene.out.tab")
        path("*.out"),                                  emit: mqc_star
    
    script:

        """
        STAR \\
            --runThreadN $task.cpus \\
            --genomeDir $params.star_genome_dir \\
            --readFilesIn $fastq_file_1 $fastq_file_2 \\
            --readFilesCommand zcat \\
            --outFilterMultimapNmax 20 \\
            --genomeLoad LoadAndRemove \\
            --outFileNamePrefix ${name}. \\
            --outSAMtype BAM SortedByCoordinate \\
            --quantMode GeneCounts \\
            --limitBAMsortRAM 20000000000
        """
}

//-----------------------------------------------------
// Post-Alignment QC (POSTQC)
//-----------------------------------------------------

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

    publishDir "$params.output/$name/alignment_post_qc/", mode: "copy", pattern: "*.{bam,bam.bai}"

    input:
        val(name)
        path(bam_file)

    output:
        val(name),                          emit: name
        path("${name}.bam"),                emit: bam_file
        path("${name}.bam.bai"),            emit: bam_index_file
        path("${name}.MarkDuplicates.txt"), emit: mqc_picard_mark_duplicates

    script:

        """
        picard -Xmx${task.memory.toGiga()}g MarkDuplicates \\
            INPUT=$bam_file \\
            OUTPUT=${name}.bam \\
            ASSUME_SORTED=true \\
            REMOVE_DUPLICATES=false \\
            METRICS_FILE=${name}.MarkDuplicates.txt \\
            VALIDATION_STRINGENCY=LENIENT \\
            TMP_DIR=./tmp/

        samtools index -@ $task.cpus ${name}.bam
        """
}

process POSTQC_STATS {

    errorStrategy "retry"
    maxRetries 3

    label "samtools"

    input:
        val(name)
        path(bam_file)
        path(bam_index_file)

    output:
        val(name),                              emit: name
        path("*.{stats,idxstats,flagstat}"),    emit: mqc_post_stats

    script:

        """
        samtools index -@ $task.cpus $bam_file
        samtools idxstats -@ $task.cpus $bam_file > ${name}.idxstats
        samtools flagstat -@ $task.cpus $bam_file > ${name}.flagstat
        samtools stats -@ $task.cpus $bam_file > ${name}.stats
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
        path("stats/*")

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
    
    A_STAR(align.name, align.fastq_1, align.fastq_2)

    //-----------------------------------------------------
    // Post-Alignment QC (POSTQC)
    //-----------------------------------------------------

    // Post-Alignment Metrics
    POSTQC_PICARD_METRICS(A_STAR.out.name_picard_metrics, A_STAR.out.bam_file_picard_metrics)

    // Mark Duplicates
    POSTQC_PICARD_MARK_DUPLICATES(A_STAR.out.name_picard_mark_duplicates, A_STAR.out.bam_file_picard_mark_duplicates)

    // Determine statistics for alignment
    POSTQC_STATS(POSTQC_PICARD_MARK_DUPLICATES.out.name, POSTQC_PICARD_MARK_DUPLICATES.out.bam_file, POSTQC_PICARD_MARK_DUPLICATES.out.bam_index_file)

    //-----------------------------------------------------
    // Summarization (SUMMARY)
    //-----------------------------------------------------

    SUMMARY_MULTI_QC(
        A_STAR.out.mqc_star.collect(),
        POSTQC_PICARD_METRICS.out.mqc_picard_metrics.collect(),
        POSTQC_PICARD_MARK_DUPLICATES.out.mqc_picard_mark_duplicates.collect(),
        POSTQC_STATS.out.mqc_post_stats.collect()
    )
}
