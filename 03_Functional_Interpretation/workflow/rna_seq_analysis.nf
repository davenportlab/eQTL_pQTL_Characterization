nextflow.enable.dsl = 2

params.samples_meta_data = "/nfs/users/nfs_n/nm18/eQTL_pQTL_Characterization/03_Functional_Interpretation/metadata/reads_calderon_et_al_rna_seq.txt"
params.rna_seq_dir = "/nfs/users/nfs_n/nm18/gains_team282/epigenetics/calderon_et_al/processed/rna_seq/"
params.rna_seq_qc_exec = "/nfs/users/nfs_n/nm18/gains_team282/nikhil/bin/rnaseqc"
params.collapsed_genome_annotation = "/nfs/users/nfs_n/nm18/eQTL_pQTL_Characterization/03_Functional_Interpretation/data/Homo_sapiens.GRCh38.99.collapsed.gtf"
params.output_dir = "/nfs/users/nfs_n/nm18/gains_team282/epigenetics/calderon_et_al/analysis/rna_seq/"

//-----------------------------------------------------
// Aggregate Fragment Counts
//-----------------------------------------------------

process AGGREGATE_STAR_FRAGMENT_COUNTS {

    errorStrategy "retry"
    maxRetries 3

    label "simple_bash"

    publishDir "$params.output_dir/", mode: "move"

    output:
        path("star_gene_counts.tsv")

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

        paste -d \$'\t' genes.tsv *.counts.tsv > star_gene_counts.tsv
        """
}

process AGGREGATE_FEATURE_COUNTS_FRAGMENT_COUNTS {

    errorStrategy "retry"
    maxRetries 3

    label "simple_bash"

    publishDir "$params.output_dir/", mode: "move"

    output:
        path("feature_counts_gene_counts.tsv")

    script:

        """
        awk -F ',' 'NR > 1 { print \$1; }' $params.samples_meta_data > samples.txt

        while read sample
        do
            echo "\$sample" > \${sample}.counts.tsv
            awk 'NR > 2 { print \$7; }' $params.rna_seq_dir/\$sample/alignment_post_qc/\${sample}.counts.tsv >> \${sample}.counts.tsv
        done <samples.txt

        first_sample=\$(head -n 1 samples.txt)

        echo "Gene" > genes.tsv
        awk 'NR > 2 { print \$1; }' $params.rna_seq_dir/\$first_sample/alignment_post_qc/\${first_sample}.counts.tsv >> genes.tsv

        paste -d \$'\t' genes.tsv *.counts.tsv > feature_counts_gene_counts.tsv
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
            $params.rna_seq_dir/$sample/alignment_post_qc/${sample}.duplicates.bam \\
            .
        
        sed -i 's/.duplicates.bam//' ${sample}.duplicates.bam.metrics.tsv
        mv ${sample}.duplicates.bam.metrics.tsv ${sample}.metrics.tsv
        """
}

process AGGREGATE_RNA_SEQC {

    errorStrategy "retry"
    maxRetries 3

    label "simple_bash"

    publishDir "$params.output_dir/", mode: "move"

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
// Workflow
//-----------------------------------------------------

workflow {

    AGGREGATE_STAR_FRAGMENT_COUNTS()

    AGGREGATE_FEATURE_COUNTS_FRAGMENT_COUNTS()

    samples_channel = Channel
        .fromPath("$params.samples_meta_data")
        .splitCsv(header: true)
        .map{ row -> row.Run }

    RNA_SEQC(samples_channel)

    AGGREGATE_RNA_SEQC(RNA_SEQC.out.rnaseqc_metrics.collect())
}
