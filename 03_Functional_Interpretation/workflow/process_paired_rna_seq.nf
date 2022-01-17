nextflow.enable.dsl = 2

params.reads_dir = "/nfs/users/nfs_n/nm18/gains_team282/epigenetics/calderon_et_al/raw/rna_seq/"
params.output_dir = "/nfs/users/nfs_n/nm18/gains_team282/epigenetics/calderon_et_al/processed/rna_seq/"
params.star_genome_dir = "/nfs/users/nfs_n/nm18/gains_team282/epigenetics/star_genome_index/"
params.genome_annotation = "/lustre/scratch118/humgen/resources/rna_seq_genomes/Homo_sapiens.GRCh38.99.gtf"
params.genome_fasta = "/lustre/scratch118/humgen/resources/rna_seq_genomes/Homo_sapiens.GRCh38.dna.primary_assembly.fa"

//-----------------------------------------------------
// Pre-Alignment QC (PREQC)
//-----------------------------------------------------

process PREQC_FASTQC {

    errorStrategy "retry"
    maxRetries 3

    label "fastqc"

    input:
        val(name)
        file(fastq_file_1)
        file(fastq_file_2)

    output:
        path("*.{zip,html}"), emit: mqc_fastqc

    script:

        def cores = (task.cpus as int) / 8

        """
        fastqc $fastq_file_1 $fastq_file_2 \\
            --threads $cores \\
            --quiet
        """
}

process PREQC_TRIMGALORE {

    errorStrategy "retry"
    maxRetries 3

    label "trim_galore"

    input:
        val(name)
        file(fastq_file_1)
        file(fastq_file_2)

    output:
        val(name),                          emit: name
        path("${name}_1.trimmed.fastq.gz"), emit: fastq_file_1
        path("${name}_2.trimmed.fastq.gz"), emit: fastq_file_2
        path("*.txt"),                      emit: mqc_trim_galore
        path("*.{zip,html}"),               emit: mqc_trim_galore_fastqc

    script:
        
        """
        trim_galore $fastq_file_1 $fastq_file_2 \\
            --cores 1 \\
            --paired \\
            --fastqc \\
            --gzip

        mv ${name}_1_val_1.fq.gz ${name}_1.trimmed.fastq.gz
        mv ${name}_2_val_2.fq.gz ${name}_2.trimmed.fastq.gz
        """
}

//-----------------------------------------------------
// Alignment (A)
//-----------------------------------------------------

process A_STAR {

    errorStrategy "retry"
    maxRetries 3

    label "star"

    publishDir "$params.output_dir/$name/alignment/", mode: "copy", pattern: "*.ReadsPerGene.out.tab"

    input:
        val(name)
        path(fastq_file_1)
        path(fastq_file_2)

    output:
        val(name),                                      emit: name_count_reads
        val(name),                                      emit: name_picard_metrics
        val(name),                                      emit: name_picard_mark_duplicates
        path("${name}.Aligned.sortedByCoord.out.bam"),  emit: bam_file_count_reads
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

process POSTQC_COUNT_READS {

    errorStrategy "retry"
    maxRetries 3

    label "featureCounts"

    publishDir "$params.output_dir/$name/alignment_post_qc/", mode: "copy"

    input:
        val(name)
        path(bam_file)
    
    output:
        path("${name}.counts.tsv"),         emit: feature_counts
        path("${name}.counts.tsv.summary"), emit: mqc_feature_counts
    
    script:

        """
        featureCounts \\
            -p \\
            -a $params.genome_annotation \\
            -t exon \\
            -g gene_id \\
            -T $task.cpus \\
            -o ${name}.counts.tsv \\
            $bam_file
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

    publishDir "$params.output_dir/$name/alignment_post_qc/", mode: "copy", pattern: "*.{bam,bam.bai}"

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
        picard -Xmx${task.memory.toGiga()}g MarkDuplicates \\
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

    publishDir "$params.output/multiqc/", mode: "move"

    input:
        path("fastqc/*")
        path("trim_galore/*")
        path("trim_galore/fastqc/*")
        path("alignment/*")
        path("feature_counts/*")
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

    qc = Channel
        .fromFilePairs("$params.reads_dir/*_{1,2}.fastq.gz")
        .multiMap { name, files ->
            name: name
            fastq_1: files[0]
            fastq_2: files[1]
        }

    align = Channel
        .fromFilePairs("$params.reads_dir/*_{1,2}.fastq.gz")
        .multiMap { name, files ->
            name: name
            fastq_1: files[0]
            fastq_2: files[1]
        }

    //-----------------------------------------------------
    // Pre-Alignment QC (PREQC)
    //-----------------------------------------------------

    // Assess Read Quality
    PREQC_FASTQC(qc.name, qc.fastq_1, qc.fastq_2)

    // Trim Low Quality Reads
    PREQC_TRIMGALORE(align.name, align.fastq_1, align.fastq_2)

    //-----------------------------------------------------
    // Alignment (A)
    //-----------------------------------------------------
    
    A_STAR(PREQC_TRIMGALORE.out.name, PREQC_TRIMGALORE.out.fastq_file_1, PREQC_TRIMGALORE.out.fastq_file_2)

    //-----------------------------------------------------
    // Post-Alignment QC (POSTQC)
    //-----------------------------------------------------

    // Count Fragments
    POSTQC_COUNT_READS(A_STAR.out.name_count_reads, A_STAR.out.bam_file_count_reads)

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
        PREQC_FASTQC.out.mqc_fastqc.collect(),
        PREQC_TRIMGALORE.out.mqc_trim_galore.collect(),
        PREQC_TRIMGALORE.out.mqc_trim_galore_fastqc.collect(),
        A_STAR.out.mqc_star.collect(),
        POSTQC_COUNT_READS.out.mqc_feature_counts.collect(),
        POSTQC_PICARD_METRICS.out.mqc_picard_metrics.collect(),
        POSTQC_PICARD_MARK_DUPLICATES.out.mqc_picard_mark_duplicates.collect(),
        POSTQC_STATS.out.mqc_post_stats.collect()
    )
}
