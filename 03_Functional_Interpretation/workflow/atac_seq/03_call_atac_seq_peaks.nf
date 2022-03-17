nextflow.enable.dsl = 2

params.atac_seq_dir = "/nfs/users/nfs_n/nm18/gains_team282/epigenetics/accessibility/processed/atac_seq/"
params.metadata = "/nfs/users/nfs_n/nm18/eQTL_pQTL_Characterization/03_Functional_Interpretation/metadata/reads_atac_seq.txt"
params.output = "/nfs/users/nfs_n/nm18/gains_team282/epigenetics/accessibility/merged/atac_seq/"


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
        path("${name}.sortedByName.bam")

    script:

        """
        samtools merge -@ $task.cpus ${name}.bam input_alignments/*.bam
        samtools index -@ $task.cpus ${name}.bam

        # Sort reads by name, making it easier to count using featureCounts
        samtools sort -m 1G -@ $task.cpus ${name}.bam > ${name}.sortedByName.bam
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
// Main Workflow
//-----------------------------------------------------

workflow {

    //-----------------------------------------------------
    // Setup
    //-----------------------------------------------------

    merge = Channel
        .fromPath(params.metadata)
        .splitCsv(skip: 1)
        .map{row -> [row[1], "$params.atac_seq_dir/${row[2]}/alignment_post_qc/${row[2]}.filtered.bam"] }
        .groupTuple()

    //-----------------------------------------------------
    // Merge Replicates (MERGE)
    //-----------------------------------------------------

    MERGE_REPLICATES(merge)

    //-----------------------------------------------------
    // Peak Calling (PEAK)
    //-----------------------------------------------------

    PEAK_MACS2_NARROW(MERGE_REPLICATES.out.name, MERGE_REPLICATES.out.merged_bam_file, MERGE_REPLICATES.out.merged_bam_index_file)
}