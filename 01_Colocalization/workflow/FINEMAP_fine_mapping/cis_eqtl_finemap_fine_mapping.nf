nextflow.enable.dsl = 2

params.output_dir = "/nfs/users/nfs_n/nm18/gains_team282/nikhil/colocalization/cis_eqtl/fine_mapping/FINEMAP/"
params.chr = "1"

process SPLIT_LOCI {

    label "Rhugemem"

    output:
        path("full/*.{master,z,ld}"),           emit: full_fm_files
        path("conditional/*.{master,z,ld}"),    emit: conditional_fm_files
    
    script:

        """
        mkdir full/
        mkdir conditional/

        Rscript $workflow.projectDir/cis_eqtl_finemap_fine_mapping_split_loci.R $params.chr
        """
}

process FINE_MAPPING_FULL {

    errorStrategy "retry"
    maxRetries 3

    label "finemap"

    input:
        path(full_fm_files)

    output:
        path("*.{config,cred,snp}*"), emit: fm_output_files

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

process FINE_MAPPING_CONDITIONAL {

    errorStrategy "retry"
    maxRetries 3

    label "finemap"

    input:
        path(conditional_fm_files)

    output:
        path("*.{config,cred,snp}*"), emit: fm_output_files

    script:

        """
        ls -1 | grep ".*\\.master" | sed s/\\.master//g > loci

        while read locus; do

            finemap \\
                --sss \\
                --in-files \${locus}.master \\
                --n-causal-snps 1 \\
                --log \\
                --n-threads $task.cpus
        
        done <loci
        """
}

process AGGREGATE_CREDIBILE_SETS {

    errorStrategy "retry"
    maxRetries 3

    publishDir "$params.output_dir/full/", mode: "move", pattern: "full_*"
    publishDir "$params.output_dir/conditional/", mode: "move", pattern: "conditional_*"

    label "R"

    input:
        path("full_cred_sets/*")
        path("conditional_cred_sets/*")

    output:
        path("full_chr${params.chr}_credible_sets.tsv")
        path("full_chr${params.chr}_pips.tsv")
        path("conditional_chr${params.chr}_credible_sets.tsv")
        path("conditional_chr${params.chr}_pips.tsv")

    script:

        """
        python3 $workflow.projectDir/cis_eqtl_finemap_fine_mapping_aggregate.py $params.chr
        """
}


workflow {

    SPLIT_LOCI()

    FINE_MAPPING_FULL(SPLIT_LOCI.out.full_fm_files)

    FINE_MAPPING_CONDITIONAL(SPLIT_LOCI.out.conditional_fm_files)

    AGGREGATE_CREDIBILE_SETS(
        FINE_MAPPING_FULL.out.fm_output_files,
        FINE_MAPPING_CONDITIONAL.out.fm_output_files
    )
}
