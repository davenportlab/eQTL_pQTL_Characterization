nextflow.enable.dsl = 2

params.reads_dir = "/nfs/users/nfs_n/nm18/gains_team282/epigenetics/accessibility/raw/rna_seq/"
params.output = "/nfs/users/nfs_n/nm18/gains_team282/epigenetics/accessibility/processed/rna_seq/"

//-----------------------------------------------------
// Pre-Alignment QC (PREQC)
//-----------------------------------------------------

process PREQC_TRIMGALORE {

    errorStrategy "retry"
    maxRetries 3

    label "trim_galore"

    publishDir "$params.output/$name/pre_processing/", mode: "move", pattern: "*.fastq.gz"

    input:
        val(name)
        file(fastq_file_1)
        file(fastq_file_2)

    output:
        path("${name}_1.trimmed.fastq.gz")
        path("${name}_2.trimmed.fastq.gz")
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
// Summarization (SUMMARY)
//-----------------------------------------------------

process SUMMARY_MULTI_QC {

    errorStrategy "retry"
    maxRetries 3

    label "multiqc"

    publishDir "$params.output/multiqc/pre_processing/", mode: "move"

    input:
        path("trim_galore/*")
        path("trim_galore/fastqc/*")

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

    //-----------------------------------------------------
    // Pre-Alignment QC (PREQC)
    //-----------------------------------------------------

    // Trim Low Quality Reads
    PREQC_TRIMGALORE(qc.name, qc.fastq_1, qc.fastq_2)
    
    //-----------------------------------------------------
    // Summarization (SUMMARY)
    //-----------------------------------------------------

    SUMMARY_MULTI_QC(
        PREQC_TRIMGALORE.out.mqc_trim_galore.collect(),
        PREQC_TRIMGALORE.out.mqc_trim_galore_fastqc.collect()
    )
}
