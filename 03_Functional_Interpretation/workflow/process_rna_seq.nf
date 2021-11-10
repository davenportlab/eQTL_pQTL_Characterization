nextflow.enable.dsl=2

params.reads_dir = "~/gains_team282/epigenetics/calderon_et_al/raw/rna_seq/"
params.output = "~/gains_team282/epigenetics/calderon_et_al/processed/rna_seq/"
params.genome_dir = "~/gains_team282/epigenetics/genome_index/"

process ReadQC {

    label "few_cpus"

    input:
        file(fastq_file)

    output:
        path("./fastqc/${base_name}_fastqc/fastqc_data.txt"), emit: fastqc_data

    script:

        base_name = fastq_file.getSimpleName()

        """
        mkdir -p ./fastqc/

        fastqc $fastq_file \
            --outdir=./fastqc/ \
            --threads 8

        unzip ./fastqc/${base_name}_fastqc.zip -d ./fastqc/

        mkdir -p $params.output/$base_name/
        cp ./fastqc/${base_name}_fastqc/fastqc_data.txt $params.output/$base_name/fastqc_data.txt
        """
}

process Trimming {

    label "few_cpus"

    input:
        file(fastq_file)

    output:
        path("./trim_galore/trimmed.fastq"), emit: fastq_file
        path("./trim_galore/${base_name}.fastq_trimming_report.txt"), emit: trim_galore_log
        val(base_name), emit: base_name

    script:

        base_name = fastq_file.getSimpleName()
        
        """
        mkdir -p ./trim_galore/

        trim_galore $fastq_file \
            --illumina \
            --output_dir ./trim_galore/ \
            --cores 8

        mv ./trim_galore/${base_name}_trimmed.fq ./trim_galore/trimmed.fastq
        """
}

process Align {

    label "big_mem"

    input:
        path(fastq_file)
        val(base_name)

    output:
        path("./star/${base_name}.Aligned.out.sam"), emit: sam_file
        path("./star/${base_name}.Log.final.out"), emit: final_log
        path("./star/${base_name}.SJ.out.tab"), emit: splicing_log
        val(base_name), emit: base_name

    script:

        """
        mkdir -p ./star/

        STAR \
            --runThreadN 12 \
            --genomeDir $params.genome_dir \
            --genomeLoad LoadAndRemove \
            --readFilesIn $fastq_file \
            --outFileNamePrefix ./star/${base_name}.

        mkdir $params.output/$base_name/
        cp ./star/${base_name}.SJ.out.tab $params.output/$base_name/${base_name}.SJ.out.tab
        """
}

process IndexAlignment {

    input:
        path(sam_file)
        val(base_name)

    output:
        path("./index/${base_name}.sorted.bam"), emit: sorted_bam
        path("./index/${base_name}.sorted.bam.idx"), emit: sorted_bam_index
        val(base_name), emit: base_name

    script:

        """
        mkdir -p ./index/

        samtools view --bam --uncompressed -@ 7 $sam_file | samtools sort -@ 7 -o ./index/${base_name}.sorted.bam -O bam
        samtools index -@ 15 ./index/${base_name}.sorted.bam ./index/${base_name}.sorted.bam.idx

        mkdir $params.output/$base_name/
        cp ./index/${base_name}.sorted.bam $params.output/$base_name/${base_name}.sorted.bam
        cp ./index/${base_name}.sorted.bam.idx $params.output/$base_name/${base_name}.sorted.bam.idx
        """
}

process MarkDuplicates {

    label "fewer_cpus_with_mem"

    input:
        path(bam_file)
        val(base_name)

    output:
        path("./duplicates/${base_name}.bam"), emit: bam_file
        path("./duplicates/mark_duplicates_metrics.txt"), emit: mark_duplicates_log
        val(base_name), emit: base_name

    script:

        """
        mkdir -p ./duplicates/

        picard MarkDuplicates \
            I=$bam_file \
            O=./duplicates/${base_name}.bam \
            M=./duplicates/mark_duplicates_metrics.txt

        mkdir $params.output/$base_name/
        cp /duplicates/${base_name}.bam $params.output/$base_name/${base_name}.duplicates_marked.bam
        """
}

process PostQC {

    label "fewer_cpus_with_mem"

    input:
        file(bam_file)
        val(base_name)
    
    output:
        path("./RSeQC/${base_name}_bam_statistics.txt"), emit: bam_statistics
        path("./RSeQC/${base_name}.GC.tsv"), emit: bam_gc
        val(base_name), emit: base_name
        
    script:

        """
        mkdir -p ./RSeQC/

        bam_stat.py -i $bam_file >> ./RSeQC/${base_name}_bam_statistics.txt

        read_GC.py -i $bam_file -o ./RSeQC/${base_name}
        mv ./RSeQC/${base_name}.GC.xls ./RSeQC/${base_name}.GC.tsv
        """
}

process CountReads {

    label "few_cpus_with_mem"

    input:
        file(bam_file)

    output:
        path("./featureCounts/counts.txt"), emit: count_matrix
        path("./featureCounts/counts.txt.summary"), emit: feature_counts_log

    script:

        """
        mkdir -p ./featureCounts/

        featureCounts \
            -T 8 \
            -s 2 \
            -t exon \
            -g gene_id \
            -a /lustre/scratch118/humgen/resources/rna_seq_genomes/Homo_sapiens.GRCh38.99.gtf \
            -o ./featureCounts/counts.txt \
            ${bam_file}

        cp ./featureCounts/counts.txt $params.output/featureCounts.tsv
        """
}

process Summarize {

    label "fewer_cpus_with_mem"

    input:
        val(base_name)
        path(fastqc_data)
        path(trim_galore_log)
        path(align_final_log)
        path(mark_duplicates_log)
        path(bam_statistics)
        path(bam_gc)
        path(feature_counts_log)

    script:

        base_names_str = base_name.join(",")

        """
        python $workflow.projectDir/process_rna_seq_summarize.py fastqc $fastqc_data $params.output/fastqc_summary.csv $base_names_str
        python $workflow.projectDir/process_rna_seq_summarize.py trim-galore $trim_galore_log $params.output/trim_galore_summary.csv $base_names_str
        python $workflow.projectDir/process_rna_seq_summarize.py star-alignment $align_final_log $params.output/align_final_summary.csv $base_names_str
        python $workflow.projectDir/process_rna_seq_summarize.py mark-duplicates $mark_duplicates_log $params.output/mark_duplicates_summary.csv $base_names_str
        python $workflow.projectDir/process_rna_seq_summarize.py bam-statistics $bam_statistics $params.output/bam_statistics_summary.csv $base_names_str
        python $workflow.projectDir/process_rna_seq_summarize.py bam-gc $bam_gc $params.output/bam_gc_summary.csv $base_names_str
        python $workflow.projectDir/process_rna_seq_summarize.py feature-counts $feature_counts_log $params.output/feature_counts_summary.csv $base_names_str
        """
}

workflow {

    fastq_channel_qc = Channel.fromPath(["$params.reads_dir/SRR7647819.fastq", "$params.reads_dir/SRR7647818.fastq"])

    ReadQC(fastq_channel_qc)

    fastq_channel_align = Channel.fromPath(["$params.reads_dir/SRR7647819.fastq", "$params.reads_dir/SRR7647818.fastq"])
    
    Trimming(fastq_channel_align)
    Align(Trimming.out.fastq_file, Trimming.out.base_name)
    IndexAlignment(Align.out.sam_file, Align.out.base_name)
    MarkDuplicates(IndexAlignment.out.sorted_bam, IndexAlignment.out.base_name)
    PostQC(MarkDuplicates.out.bam_file, MarkDuplicates.out.base_name)

    CountReads(MarkDuplicates.out.bam_file.collect())

    Summarize(
        PostQC.out.base_name.collect(),
        ReadQC.out.fastqc_data.collect(),
        Trimming.out.trim_galore_log.collect(),
        Align.out.final_log.collect(),
        MarkDuplicates.out.mark_duplicates_log.collect(),
        PostQC.out.bam_statistics.collect(),
        PostQC.out.bam_gc.collect(),
        CountReads.out.feature_counts_log.collect()
    )
}