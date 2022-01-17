nextflow.enable.dsl = 2

params.samples_meta_data = "/nfs/users/nfs_n/nm18/eQTL_pQTL_Characterization/03_Functional_Interpretation/metadata/reads_calderon_et_al_rna_seq.txt"
params.rna_seq_dir = "/nfs/users/nfs_n/nm18/gains_team282/epigenetics/calderon_et_al/processed/rna_seq/"
params.genome_annotation = "/lustre/scratch118/humgen/resources/rna_seq_genomes/Homo_sapiens.GRCh38.99.gtf"
params.output_dir = "/nfs/users/nfs_n/nm18/gains_team282/epigenetics/calderon_et_al/analysis/rna_seq/"

//-----------------------------------------------------
// Aggregate Fragment Counts
//-----------------------------------------------------

process AGGREGATE_FRAGMENT_COUNTS {

    errorStrategy "retry"
    maxRetries 3

    label "simple_bash"

    publishDir "$params.output_dir/", mode: "move"

    output:
        path("gene_counts.tsv")

    script:

        """
        awk -F ',' 'NR > 1 { print \$1; }' $params.samples_meta_data > samples.txt

        while read sample
        do
            echo "\$sample" > \${sample}.counts.tsv
            awk 'NR > 4 { print \$2; }' $params.rna_seq_dir/\$sample/alignment/\${sample}.ReadsPerGene.out.tab >> \${sample}.counts.tsv
        done <samples.txt

        first_sample=\$(head -n 1 samples.txt)

        echo "Gene" > genes.tsv
        awk 'NR > 4 { print \$1; }' $params.rna_seq_dir/\$first_sample/alignment/\${first_sample}.ReadsPerGene.out.tab >> genes.tsv

        paste -d \$'\t' genes.tsv *.counts.tsv > gene_counts.tsv
        """
}

//-----------------------------------------------------
// RNA-Seq QC
//-----------------------------------------------------

process RNA_SEQ_QUALIMAP {

    errorStrategy "retry"
    maxRetries 3

    label "qualimap"

    input:
        val(sample)

    output:
        path("bamqc/data.tsv"),     emit: bam_qc_data
        path("rnaseq/data.tsv"),    emit: rna_seq_qc_data
    
    script:

        """
        qualimap bamqc \\
            -bam $params.rna_seq_dir/$sample/alignment_post_qc/${sample}.duplicates.bam \\
            -outdir ./bamqc/ \\
            -nt $task.cpus \\
            --java-mem-size=10G

        qualimap rnaseq \\
            -bam $params.rna_seq_dir/$sample/alignment_post_qc/${sample}.duplicates.bam \\
            -gtf $params.genome_annotation \\
            --paired --sorted \\
            -outdir ./rnaseq/ \\
            --java-mem-size=10G

        # Extract data from BAM QC

        header=""
        qc_data=""
        
        while read -r line
        do
            key=\$(echo "\$line" | cut -d "=" -f 1 | awk '{\$1=\$1};1')
            value=\$(echo "\$line" | cut -d "=" -f 2 | awk '{\$1=\$1};1')
            header="\${header}\t\\"\$key\\""
            qc_data="\${qc_data}\t\\"\$value\\""
        done <<<\$(grep "=" ./bamqc/genome_results.txt)

        echo -e "\$header" > ./bamqc/data.tsv
        echo -e "\$qc_data" >> ./bamqc/data.tsv

        # Extract data from RNA-Seq QC

        header=""
        qc_data=""
        
        while read -r line
        do
            key=\$(echo "\$line" | cut -d "=" -f 1 | awk '{\$1=\$1};1')
            value=\$(echo "\$line" | cut -d "=" -f 2 | awk '{\$1=\$1};1')
            header="\${header}\t\\"\$key\\""
            qc_data="\${qc_data}\t\\"\$value\\""
        done <<<\$(grep "=" ./rnaseq/rnaseq_qc_results.txt)

        echo -e "\$header" > ./rnaseq/data.tsv
        echo -e "\$qc_data" >> ./rnaseq/data.tsv
        """
}

process AGGREGATE_RNA_SEQ_QC {

    errorStrategy "retry"
    maxRetries 3

    label "simple_bash"

    publishDir "$params.output_dir/", mode: "move"

    input:
        path("bamqc/*.tsv")
        path("rnaseq/*.tsv")

    output:
        path("bam_qc.tsv")
        path("rna_seq_qc.tsv")
    
    script:

        """
        awk -F ',' 'NR > 1 { print \$1; }' $params.samples_meta_data > samples.txt

        # BAM QC Data File

        first_file=\$(ls -1 bamqc/*.tsv | head -n 1)
        header_line=\$(head -n 1 \$first_file)
        echo -e "Sample\$header_line" > bam_qc.tsv

        awk 'FNR>1' bamqc/*.tsv > bam_qc_data.tsv

        paste -d '' samples.txt bam_qc_data.tsv >> bam_qc.tsv

        # RNA-Seq QC Data File

        first_file=\$(ls -1 rnaseq/*.tsv | head -n 1)
        header_line=\$(head -n 1 \$first_file)
        echo -e "Sample\$header_line" > rna_seq_qc.tsv

        awk 'FNR>1' rnaseq/*.tsv > rna_seq_qc_data.tsv

        paste -d '' samples.txt rna_seq_qc_data.tsv >> rna_seq_qc.tsv
        """
}

//-----------------------------------------------------
// Workflow
//-----------------------------------------------------

workflow {

    AGGREGATE_FRAGMENT_COUNTS()

    samples_channel = Channel
        .fromPath("$params.samples_meta_data")
        .splitCsv(header: true)
        .map{ row -> row.Run }

    RNA_SEQ_QUALIMAP(samples_channel)

    AGGREGATE_RNA_SEQ_QC(
        RNA_SEQ_QUALIMAP.out.bam_qc_data.collect(),
        RNA_SEQ_QUALIMAP.out.rna_seq_qc_data.collect()
    )
}
