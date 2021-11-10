nextflow.enable.dsl = 2

params.reads_dir = "~/gains_team282/epigenetics/calderon_et_al/raw/atac_seq/"
params.output = "~/gains_team282/epigenetics/calderon_et_al/processed/atac_seq/"
params.star_genome_dir = "~/gains_team282/epigenetics/star_genome_index/"
params.bowtie2_genome_dir = "~/gains_team282/epigenetics/bowtie2_genome_index/"
params.genome_fasta = "/lustre/scratch118/humgen/resources/rna_seq_genomes/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
params.genome_black_list = "~/eQTL_pQTL_Characterization/03_Functional_Interpretation/data/hg38-blacklist.v2.bed"

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

        """
        fastqc $fastq_file_1 $fastq_file_2 \\
            --threads $task.cpus \\
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

        def cores = (task.cpus as int) - 4
        
        """
        trim_galore $fastq_file_1 $fastq_file_2 \\
            --cores $cores \\
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

    input:
        val(name)
        path(fastq_file_1)
        path(fastq_file_2)

    output:
        val(name),                                          emit: name
        path("${name}.Aligned.sortedByCoord.out.bam"),      emit: bam_file
        path("${name}.Aligned.sortedByCoord.out.bam.bai"),  emit: bam_index_file
        path("${name}.Log.final.out"),                      emit: mqc_alignment

    script:

        """
        STAR \
            --runThreadN $task.cpus \\
            --genomeDir $params.star_genome_dir \\
            --genomeLoad LoadAndRemove \\
            --limitBAMsortRAM 20000000000 \\
            --readFilesIn $fastq_file_1 $fastq_file_2 \\
            --readFilesCommand zcat \\
            --outWigType bedGraph \\
            --outSAMtype BAM SortedByCoordinate \\
            --outSAMunmapped Within \\
            --quantMode GeneCounts \\
            --outFileNamePrefix ${name}.

        samtools index -@ $task.cpus ${name}.Aligned.sortedByCoord.out.bam
        """
}

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
        path("${name}.filtered.sorted.bam"),        emit: bam_file
        path("${name}.filtered.sorted.bam.bai"),    emit: bam_index_file
        path("*.{stats,idxstats,flagstat}"),        emit: mqc_post_filtering

    script:

        // -f 0x0001 - Keep reads that are paired.
        // -F 0x0004 - Filter out reads that are unmapped.
        // -F 0x0200 -F 0x0400 - Filter out reads that are marked as duplicates.
        // -q 1 - Only keeps reads with MAPQ Score >= 1. MAPQ is only 0 when read is multi-mapped.

        """
        samtools view $bam_file \\
            -f 0x0001 \\
            -F 0x0004 \\
            -F 0x0200 -F 0x0400 \\
            -q 1 \\
            --bam \\
            -@ $task.cpus > ${name}.filtered.bam
        
        samtools index -@ $task.cpus ${name}.filtered.bam

        # Filter reads that map to mitochondrial genome (MT)
        samtools idxstats ${name}.filtered.bam | cut -f 1 | grep -v MT | xargs samtools view -b ${name}.filtered.bam > ${name}.filtered.noMT.bam

        # Filter reads in blacklisted region from ENCODE
        bedtools intersect -abam ${name}.filtered.noMT.bam -b $params.genome_black_list -v > ${name}.filtered.bam

        samtools sort -m 1G -@ $task.cpus ${name}.filtered.bam > ${name}.filtered.sorted.bam
        samtools index -@ $task.cpus ${name}.filtered.sorted.bam
        samtools idxstats -@ $task.cpus ${name}.filtered.sorted.bam > ${name}.filtered.sorted.bam.idxstats
        samtools flagstat -@ $task.cpus ${name}.filtered.sorted.bam > ${name}.filtered.sorted.bam.flagstat
        samtools stats -@ $task.cpus ${name}.filtered.sorted.bam > ${name}.filtered.sorted.bam.stats
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

    publishDir "$params.output/multiqc/", mode: "move"

    input:
        path("fastqc/*")
        path("trim_galore/*")
        path("trim_galore/fastqc/*")
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
    
    // Option 1: STAR
    // A_STAR(PREQC_TRIMGALORE.out.name, PREQC_TRIMGALORE.out.fastq_file_1, PREQC_TRIMGALORE.out.fastq_file_2)

    // Option 2: Bowtie2
    A_BOWTIE2(PREQC_TRIMGALORE.out.name, PREQC_TRIMGALORE.out.fastq_file_1, PREQC_TRIMGALORE.out.fastq_file_2)

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
        PREQC_FASTQC.out.mqc_fastqc.collect(),
        PREQC_TRIMGALORE.out.mqc_trim_galore.collect(),
        PREQC_TRIMGALORE.out.mqc_trim_galore_fastqc.collect(),
        POSTQC_SORT.out.mqc_alignment.collect(),
        POSTQC_PICARD_METRICS.out.mqc_picard_metrics.collect(),
        POSTQC_PICARD_MARK_DUPLICATES.out.mqc_picard_mark_duplicates.collect(),
        POSTQC_FILTERING.out.mqc_post_filtering.collect()
    )
}
