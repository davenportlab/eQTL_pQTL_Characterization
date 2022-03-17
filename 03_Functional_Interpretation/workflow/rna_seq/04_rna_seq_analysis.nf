nextflow.enable.dsl = 2

params.metadata = "/nfs/users/nfs_n/nm18/eQTL_pQTL_Characterization/03_Functional_Interpretation/metadata/reads_rna_seq.txt"
params.rna_seq_dir = "/nfs/users/nfs_n/nm18/gains_team282/epigenetics/accessibility/merged/rna_seq/"
params.rna_seq_qc_exec = "/nfs/users/nfs_n/nm18/gains_team282/nikhil/bin/rnaseqc"
params.collapsed_genome_annotation = "/nfs/users/nfs_n/nm18/eQTL_pQTL_Characterization/03_Functional_Interpretation/data/Homo_sapiens.GRCh38.99.collapsed.gtf"
params.output = "/nfs/users/nfs_n/nm18/gains_team282/epigenetics/accessibility/analysis/rna_seq/"

//-----------------------------------------------------
// Aggregate Fragment Counts
//-----------------------------------------------------

process AGGREGATE_GENE_COUNTS {

    errorStrategy "retry"
    maxRetries 3

    label "simple_bash"

    publishDir "$params.output/", mode: "move"

    output:
        path("gene_counts.tsv")

    script:

        """
        awk -F ',' 'NR > 1 { print \$2; }' $params.metadata | sort | uniq > samples.txt

        while read sample
        do
            echo "\$sample" > \${sample}.counts.tsv
            awk 'NR > 2 { print \$7; }' $params.rna_seq_dir/\$sample/counts/\${sample}.counts.tsv >> \${sample}.counts.tsv
        done <samples.txt

        first_sample=\$(head -n 1 samples.txt)

        echo "Gene" > genes.tsv
        awk 'NR > 2 { print \$1; }' $params.rna_seq_dir/\$first_sample/counts/\${first_sample}.counts.tsv >> genes.tsv

        paste -d \$'\t' genes.tsv *.counts.tsv > gene_counts.tsv
        """
}

//-----------------------------------------------------
// RNA-Seq QC
//-----------------------------------------------------

process RNA_SEQC {

    errorStrategy "retry"
    maxRetries 3

    label "rnaseqc"

    input:
        val(sample)
    
    output:
        path("${sample}.metrics.tsv"), emit: rnaseqc_metrics
        
    script:

        """
        $params.rna_seq_qc_exec \\
            $params.collapsed_genome_annotation \\
            $params.rna_seq_dir/$sample/alignment/${sample}.bam \\
            .
        
        sed -i 's/.bam//' ${sample}.bam.metrics.tsv
        mv ${sample}.bam.metrics.tsv ${sample}.metrics.tsv
        """
}

process AGGREGATE_RNA_SEQC {

    errorStrategy "retry"
    maxRetries 3

    label "simple_bash"

    publishDir "$params.output/", mode: "move"

    input:
        path("qc/*.tsv")
    
    output:
        path("rnaseqc.tsv")
    
    script:
        
        """
        mkdir qc_clean/

        for file in qc/*.tsv
        do
            file=\$(basename \$file)
            awk -F \$'\t' '{ print \$2; }' qc/\$file > qc_clean/\$file
        done

        first_sample=\$(ls qc/ | head -n 1)

        echo "Metric" > metrics.tsv
        awk -F \$'\t' 'NR > 1 { print \$1; }' qc/\$first_sample >> metrics.tsv

        paste -d \$'\t' metrics.tsv qc_clean/*.tsv > rnaseqc.tsv
        """
}


//-----------------------------------------------------
// Main Workflow
//-----------------------------------------------------

workflow {

    AGGREGATE_GENE_COUNTS()

    samples_channel = Channel
        .fromPath(params.metadata)
        .splitCsv(skip: 1)
        .map{row -> row[1]}

    RNA_SEQC(samples_channel)

    AGGREGATE_RNA_SEQC(RNA_SEQC.out.rnaseqc_metrics.collect())
}
