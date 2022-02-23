nextflow.enable.dsl = 2

params.n_proteins = "269"
params.pqtl_summary_dir = "/nfs/users/nfs_n/nm18/gains_team282/proteomics/pqtl/pqtl_ms2019/results/"
params.output_dir = "/nfs/users/nfs_n/nm18/gains_team282/proteomics/pqtl/pqtl_ms2019/"

process FIND_PQTL {

    label "Rbigmem"

    input:
        tuple(val(protein_name), path("summary_stats/*"))

    output:
        path("pqtl.summary.RDS"),       emit: wg_pqtl_summary
        path("cis.pqtl.summary.RDS"),   emit: cis_pqtl_summary

    script:

        """
        Rscript $workflow.projectDir/extract_pqtl.R $protein_name $params.n_proteins
        """
}

process AGGREGATE_PQTL {

    label "Rbigmem"

    publishDir "$params.output_dir/", mode: "move"

    input:
        path("wg_pqtl/*.RDS")
        path("cis_pqtl/*.RDS")
    
    output:
        path("whole_genome_pqtl_all.RDS")
        path("cis_pqtl_all.RDS")

    script:

        """
        Rscript $workflow.projectDir/extract_pqtl_aggregate.R
        """
}

workflow {

    summary_files = Channel.fromFilePairs("$params.pqtl_summary_dir/Protein_*_group_{1,2,3,4}.rds", size: -1)

    FIND_PQTL(summary_files)

    AGGREGATE_PQTL(
        FIND_PQTL.out.wg_pqtl_summary.collect(),
        FIND_PQTL.out.cis_pqtl_summary.collect()
    )
}
