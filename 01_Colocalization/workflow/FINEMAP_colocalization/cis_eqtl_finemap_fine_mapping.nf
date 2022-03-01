nextflow.enable.dsl = 2

params.output_dir = "/nfs/users/nfs_n/nm18/gains_team282/nikhil/colocalization/cis_eqtl/fine_mapping/FINEMAP/"
params.chr = "1"

process SPLIT_LOCI {

    label "Rbigmem"

    output:
        path("*.master"),   emit: master_files
        path("*.z"),        emit: z_files
        path("*.ld"),       emit: ld_files
    
    script:

        """
        Rscript $workflow.projectDir/cis_eqtl_finemap_fine_mapping_split_loci.R $params.chr
        """
}

process FINE_MAPPING {

    errorStrategy "retry"
    maxRetries 3

    label "finemap"

    input:
        path(master_files)
        path(z_files)
        path(ld_files)

    output:
        path("*.config"), emit: config_files
        path("*.cred*"), emit: credible_set_files
        path("*.snp"), emit: snp_files

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
        path("credible_sets/*")
        path("snps/*")

    output:
        path("chr${params.chr}_credible_sets.tsv")
        path("chr${params.chr}_pips.tsv")

    script:

        """
        python3 $workflow.projectDir/cis_eqtl_finemap_fine_mapping_aggregate.py $params.chr
        """
}


workflow {

    SPLIT_LOCI()

    FINE_MAPPING(
        SPLIT_LOCI.out.master_files,
        SPLIT_LOCI.out.z_files,
        SPLIT_LOCI.out.ld_files
    )

    AGGREGATE_CREDIBILE_SETS(
        FINE_MAPPING.out.credible_set_files,
        FINE_MAPPING.out.snp_files
    )
}
