nextflow.enable.dsl = 2

params.output_dir = "/nfs/users/nfs_n/nm18/gains_team282/nikhil/colocalization/pqtl/fine_mapping/FINEMAP/"

process SPLIT_LOCI {

    label "Rbigmem"

    output:
        path("cis_pqtl/*.master"),      emit: cis_master_files
        path("cis_pqtl/*.z"),           emit: cis_z_files
        path("cis_pqtl/*.ld"),          emit: cis_ld_files
        path("trans_pqtl/*.master"),    emit: trans_master_files
        path("trans_pqtl/*.z"),         emit: trans_z_files
        path("trans_pqtl/*.ld"),        emit: trans_ld_files
    
    script:

        """
        Rscript $workflow.projectDir/pqtl_finemap_fine_mapping_split_loci.R
        """
}

process FINE_MAPPING {

    errorStrategy "retry"
    maxRetries 3

    label "finemap"

    input:
        val(cis_trans)
        path(master_files)
        path(z_files)
        path(ld_files)

    output:
        val(cis_trans),     emit: cis_trans
        path("*.config"),   emit: config_files
        path("*.cred*"),    emit: credible_set_files
        path("*.snp"),      emit: snp_files

    script:

        """
        ls -1 | grep ".*\\.master" | sed s/\\.master//g > genes

        while read gene; do

            N_SNPS=\$(echo "\$(wc -l < \${gene}.z) - 1" | bc)

            N_CAUSAL_SNPS=10
            if [ \$N_SNPS -lt 10 ]; then
                N_CAUSAL_SNPS=\$N_SNPS
            fi

            finemap \\
                --sss \\
                --in-files \${gene}.master \\
                --n-causal-snps \$N_CAUSAL_SNPS \\
                --log \\
                --n-threads $task.cpus
        
        done <genes
        """
}

process AGGREGATE_CREDIBILE_SETS {

    errorStrategy "retry"
    maxRetries 3

    publishDir "$params.output_dir/", mode: "move"

    label "R"

    input:
        val(cis_trans)
        path("credible_sets/*")
        path("snps/*")

    output:
        path("*_credible_sets.tsv")
        path("*_pips.tsv")

    script:

        """
        python3 $workflow.projectDir/pqtl_finemap_fine_mapping_aggregate.py $cis_trans
        """
}


workflow {

    SPLIT_LOCI()

    cis_trans = Channel.from("cis", "trans")
    master_files = SPLIT_LOCI.out.cis_master_files.concat(SPLIT_LOCI.out.trans_master_files)
    z_files = SPLIT_LOCI.out.cis_z_files.concat(SPLIT_LOCI.out.trans_z_files)
    ld_files = SPLIT_LOCI.out.cis_ld_files.concat(SPLIT_LOCI.out.trans_ld_files)

    FINE_MAPPING(cis_trans, master_files, z_files, ld_files)

    AGGREGATE_CREDIBILE_SETS(
        FINE_MAPPING.out.cis_trans,
        FINE_MAPPING.out.credible_set_files,
        FINE_MAPPING.out.snp_files
    )
}
