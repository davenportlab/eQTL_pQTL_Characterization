nextflow.enable.dsl=2

params.section = "accessibility"
params.assay = "rna_seq"
params.sra_table = "/nfs/users/nfs_n/nm18/eQTL_pQTL_Characterization/03_Functional_Interpretation/metadata/reads_rna_seq.txt"
params.output_dir = "/nfs/users/nfs_n/nm18/gains_team282/epigenetics/"


process DUMP_FASTQ {

    errorStrategy "retry"
    maxRetries 5

    label "dump_fastq"

    publishDir "$params.output_dir/$params.section/raw/$params.assay/", mode: "move"

    input:
        val(run)
        val(donor)
        val(cell_type)
        val(treatment)
        val(replicate)
        val(layout)

    output:
        path("*.fastq.gz")

    script:

        if (layout == "SINGLE")

            """
            fasterq-dump --threads $task.cpus $run

            gzip ${run}.fastq

            mv ${run}.fastq.gz ${donor}-${cell_type}-${treatment}-${replicate}.fastq.gz
            """
        
        else if (layout == "PAIRED")

            """
            fasterq-dump --threads $task.cpus $run

            gzip ${run}_1.fastq
            gzip ${run}_2.fastq

            mv ${run}_1.fastq.gz ${donor}-${cell_type}-${treatment}-${replicate}_1.fastq.gz
            mv ${run}_2.fastq.gz ${donor}-${cell_type}-${treatment}-${replicate}_2.fastq.gz
            """
}

//-----------------------------------------------------
// Workflow
//-----------------------------------------------------

workflow {

    accessions = Channel
        .fromPath(params.sra_table)
        .splitCsv(header:true)
        .multiMap{row -> 
            run: row.Run
            donor: row.Donor
            cell_type: row.Cell_type
            treatment: row.Treatment
            replicate: row.Replicate
            layout: row.LibraryLayout
        }

    DUMP_FASTQ(
        accessions.run,
        accessions.donor,
        accessions.cell_type,
        accessions.treatment,
        accessions.replicate,
        accessions.layout
    )
}
