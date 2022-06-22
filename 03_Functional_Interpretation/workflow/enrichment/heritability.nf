nextflow.enable.dsl = 2

params.genotypes_prefix = "/nfs/users/nfs_n/nm18/gains_team282/Genotyping/All_genotyping_merged_filtered_b38_refiltered_rsID"
params.mapping_patients = "/nfs/users/nfs_n/nm18/gains_team282/nikhil/expression/eigengene_sva/mapping_patients.txt"
params.output_dir = "/nfs/users/nfs_n/nm18/gains_team282/epigenetics/enrichment/heritability/"


process GRM_CALCULATION {

    script:

        """
        plink \\
            --bfile $params.genotypes_prefix \\
            --keep $params.mapping_patients \\
            --make-bed \\
            --out ./autosomes \\
            --allow-extra-chr \\
            --maf 0.01 \\
            --chr 1-22

        gcta64 \\
            --bfile autosomes \\
            --make-grm \\
            --out full_grm \\
            --thread-num 16

        awk 'OFS="\\t" { print \$1, \$4 - 1, \$4, \$2; }' autosomes.bim > snps.bed
        

        bedtools intersect -a snps.bed -b /nfs/users/nfs_n/nm18/gains_team282/epigenetics/accessibility/analysis/atac_seq/immune/peak_sets/Monocytes-Stimulated.peaks.bed -wa | awk '{ print \$4; }' > monocytes.snps
        bedtools intersect -v -a snps.bed -b /nfs/users/nfs_n/nm18/gains_team282/epigenetics/accessibility/analysis/atac_seq/immune/peak_sets/Monocytes-Stimulated.peaks.bed -wa | awk '{ print \$4; }' > monocytes.not.snps

        gcta64 \\
            --bfile autosomes \\
            --extract monocytes.snps \\
            --make-grm \\
            --out monocytes_grm \\
            --thread-num 16

        gcta64 \\
            --bfile autosomes \\
            --extract monocytes.not.snps \\
            --make-grm \\
            --out monocytes_not_grm \\
            --thread-num 16
        """
}

workflow {

    GRM_CALCULATION()
}
