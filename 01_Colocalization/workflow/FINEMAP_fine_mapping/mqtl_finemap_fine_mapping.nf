nextflow.enable.dsl = 2

params.output_dir = "/nfs/users/nfs_n/nm18/gains_team282/nikhil/colocalization/mqtl/fine_mapping/FINEMAP/"

process SPLIT_LOCI {

    label "Rhugemem"

    output:
        path("*.{master,z,ld}"), emit: mqtl_files
    
    script:

        """
        Rscript $workflow.projectDir/mqtl_finemap_fine_mapping_split_loci.R
        """
}

process FINE_MAPPING {

    errorStrategy "retry"
    maxRetries 3

    label "finemap"

    input:
        path(mqtl_files)

    output:
        path("*.config"),   emit: config_files
        path("*.cred*"),    emit: credible_set_files
        path("*.snp"),      emit: snp_files

    script:

        """
        ls -1 | grep ".*\\.master" | sed s/\\.master//g > mqtl_loci

        while read mqtl_locus; do

            N_SNPS=\$(echo "\$(wc -l < \${mqtl_locus}.z) - 1" | bc)

            N_CAUSAL_SNPS=10
            if [ \$N_SNPS -lt 10 ]; then
                N_CAUSAL_SNPS=\$N_SNPS
            fi

            finemap \\
                --sss \\
                --in-files \${mqtl_locus}.master \\
                --n-causal-snps \$N_CAUSAL_SNPS \\
                --log \\
                --n-threads $task.cpus
        
        done <mqtl_loci
        """
}

process AGGREGATE_CREDIBILE_SETS {

    errorStrategy "retry"
    maxRetries 3

    publishDir "$params.output_dir/", mode: "move"

    label "R"

    input:
        path("credible_sets/*")
        path("snps/*")

    output:
        path("*_credible_sets.tsv")
        path("*_pips.tsv")

    script:

        """
        python3 $workflow.projectDir/mqtl_finemap_fine_mapping_aggregate.py
        """
}


workflow {

    SPLIT_LOCI()

    FINE_MAPPING(SPLIT_LOCI.out.mqtl_files)

    AGGREGATE_CREDIBILE_SETS(
        FINE_MAPPING.out.credible_set_files,
        FINE_MAPPING.out.snp_files
    )
}
