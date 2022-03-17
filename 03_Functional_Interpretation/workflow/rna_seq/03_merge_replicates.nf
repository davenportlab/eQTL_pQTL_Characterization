nextflow.enable.dsl = 2

params.rna_seq_dir = "/nfs/users/nfs_n/nm18/gains_team282/epigenetics/accessibility/processed/rna_seq/"
params.metadata = "/nfs/users/nfs_n/nm18/eQTL_pQTL_Characterization/03_Functional_Interpretation/metadata/reads_rna_seq.txt"
params.genome_annotation = "/lustre/scratch118/humgen/resources/rna_seq_genomes/Homo_sapiens.GRCh38.99.gtf"
params.output = "/nfs/users/nfs_n/nm18/gains_team282/epigenetics/accessibility/merged/rna_seq/"


//-----------------------------------------------------
// Merge Replicates (MERGE)
//-----------------------------------------------------

process MERGE_REPLICATES {

    errorStrategy "retry"
    maxRetries 3

    label "samtools"

    publishDir "$params.output/$name/alignment/", mode: "copy"
    
    input:
        tuple val(name), path("input_alignments/*")
    
    output:
        val(name),                          emit: name
        path("${name}.bam"),                emit: merged_bam_file
        path("${name}.bam.bai"),            emit: merged_bam_index_file
        path("${name}.sortedByName.bam"),   emit: merged_bam_file_sorted

    script:

        """
        samtools merge -@ $task.cpus ${name}.bam input_alignments/*.bam
        samtools index -@ $task.cpus ${name}.bam

        # Sort reads by name, making it easier to count using featureCounts
        samtools sort -m 1G -@ $task.cpus ${name}.bam > ${name}.sortedByName.bam
        """
}

//-----------------------------------------------------
// Read Counting (COUNT)
//-----------------------------------------------------

process COUNT_READS {

    errorStrategy "retry"
    maxRetries 3

    label "featureCounts"

    publishDir "$params.output/$name/counts/", mode: "move"

    input:
        val(name)
        path(bam_file)
    
    output:
        path("${name}.counts.tsv")
    
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


//-----------------------------------------------------
// Main Workflow
//-----------------------------------------------------

workflow {

    //-----------------------------------------------------
    // Setup
    //-----------------------------------------------------

    merge = Channel
        .fromPath(params.metadata)
        .splitCsv(skip: 1)
        .map{row -> [row[1], "$params.rna_seq_dir/${row[2]}/alignment_post_qc/${row[2]}.bam"] }
        .groupTuple()

    //-----------------------------------------------------
    // Merge Replicates (MERGE)
    //-----------------------------------------------------

    MERGE_REPLICATES(merge)

    //-----------------------------------------------------
    // Read Counting (COUNT)
    //-----------------------------------------------------

    COUNT_READS(MERGE_REPLICATES.out.name, MERGE_REPLICATES.out.merged_bam_file_sorted)
}
