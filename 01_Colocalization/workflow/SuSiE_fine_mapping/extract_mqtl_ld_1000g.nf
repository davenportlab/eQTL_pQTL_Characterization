nextflow.enable.dsl = 2

params.mqtl_sum_stats = "/nfs/users/nfs_n/nm18/gains_team282/nikhil/expression/eigengene_sva/wgcna_summary_statistics/"
params.vcf_dir = "/lustre/scratch118/humgen/resources/1000g/release/20130502/"
params.eur_samples = "/nfs/users/nfs_n/nm18/eQTL_pQTL_Characterization/01_Colocalization/data/1000G/EUR.samples.txt"
params.output_dir = "/nfs/users/nfs_n/nm18/gains_team282/nikhil/data/1000G/"


process GET_LOCUS_INFO {

    label "simple_bash"

    output:
        path("*.tsv"), emit: loci

    script:

        """
        ls -1 $params.mqtl_sum_stats | \\
            grep '^ME_[[:digit:]]\\+_' | \\
            sed 's/ME_[[:digit:]]\\+_[[:digit:]]-//g' | \\
            sed 's/\\.tsv//g' | sort | uniq > loci.txt

        while read locus; do

            one_mqtl_file=$params.mqtl_sum_stats/\$(ls -1 $params.mqtl_sum_stats | grep \$locus | head -n 1)
            awk '{ print \$1; }' \$one_mqtl_file > \${locus}.tsv

        done <loci.txt
        """
}

process LOCUS_LD {

    label "simple_bash"

    publishDir "$params.output_dir/", mode: "move"

    input:
        path(snps)

    output:
        path("${snps.getSimpleName()}.bim")
        path("${snps.getSimpleName()}.raw")

    script:

        """
        chr=\$(echo "${snps.getSimpleName()}" | sed 's/:.*\$//g')
        file=ALL.chr\${chr}.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz

        vcftools \\
            --gzvcf $params.vcf_dir/\$file \\
            --snps $snps \\
            --min-alleles 2 \\
            --max-alleles 2 \\
            --keep $params.eur_samples \\
            --plink \\
            --out plink_${snps.getSimpleName()}

        plink \\
            --file plink_${snps.getSimpleName()} \\
            --recode A \\
            --make-bed \\
            --out ${snps.getSimpleName()}
        """
}


workflow {

    GET_LOCUS_INFO()

    LOCUS_LD(GET_LOCUS_INFO.out.loci.flatten())
}
