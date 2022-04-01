nextflow.enable.dsl = 2

params.genotype_info = "/nfs/users/nfs_n/nm18/gains_team282/Genotyping/All_genotyping_merged_filtered_b38_refiltered_rsID.bim"
params.gene_info = "/nfs/team282/data/gains_team282/gene_info_864_20412_hla.txt"
params.lead_snps = "/nfs/users/nfs_n/nm18/gains_team282/nikhil/colocalization/cis_eqtl/fine_mapping/LD/lead_snps.80r2.tags.tsv"
params.conditional_snps = "/nfs/users/nfs_n/nm18/gains_team282/nikhil/colocalization/cis_eqtl/fine_mapping/LD/conditional_snps.80r2.tags.tsv"
params.conditional_results = "/nfs/users/nfs_n/nm18/gains_team282/eqtl/cisresults/conditionalanalysis/conditional_eQTL_results_final.txt"
params.sepsis_snps = "/lustre/scratch119/humgen/projects/gains_team282/eqtl/gtex/sepsis_specific_gains_lead_mashr.txt"
params.mashr_results = "/lustre/scratch119/humgen/projects/gains_team282/eqtl/gtex/gains_gtex_mashr_results.txt"
params.da_peaks = "/nfs/users/nfs_n/nm18/gains_team282/epigenetics/accessibility/analysis/atac_seq/immune/da_peak_set.csv"
params.output_dir = "/nfs/users/nfs_n/nm18/gains_team282/epigenetics/enrichment/go_shifter/"


process PREPARE_SNP_LIST {

    label "simple_bash"

    output:
        path("snps_lists/*.txt"), emit: snps_lists
    
    script:

        """
        mkdir snps_lists/

        ### LEAD SNPS

        awk -F "\\t" 'NR > 1 { print \$1; }' $params.lead_snps | sort | uniq > snps_lists/lead_snps.txt
        
        ### LEAD SNPS FROM CONDITIONAL ANALYSIS

        awk -F "\\t" 'NR > 1 { print \$1; }' $params.conditional_snps | sort | uniq > snps_lists/conditional_snps.txt
        
        ### LEAD SNPS FROM CONDITIONAL ANALYSIS W/O PRIMARY EFFECT

        sed 's/"//g' $params.conditional_results | awk 'NR > 1 { if (\$7 > 1) { print \$2; } }' | sort | uniq > snps_lists/conditional_secondary_snps.txt
        
        ### LEAD SNPS W/O SEPSIS-SPECIFIC SNPS

        # Extract lead SNPs identified from initial cis-eQTL pass
        awk -F "\\t" 'NR > 1 { print \$1; }' $params.lead_snps > lead_snps.txt

        # Extract SNPs from mashR results
        awk 'NR > 1 { print \$16; }' $params.sepsis_snps > sepsis_snps.txt

        # Filter lead SNPs to non-sepsis-specific eSNPs
        grep -vwFf sepsis_snps.txt lead_snps.txt | sort | uniq > snps_lists/non_sepsis_snps.txt

        ### SEPSIS-SPECIFIC SNPS

        # Extract SNPs from mashR results
        awk 'NR > 1 { print \$16; }' $params.sepsis_snps | sort | uniq > snps_lists/sepsis_snps.txt

        ### SEPSIS-SPECIFIC SNPS STRONGER IN GAinS

        # Get SNP regions that are stronger in sepsis from the mashR results
        awk 'NR > 1 { if (\$5 == "TRUE" && (\$2)^2 > (\$3)^2) { print \$1; } }' $params.mashr_results | sort > sepsis_up_snp_regions.txt

        # Filter sepsis snps based on mashR results
        grep -wFf sepsis_up_snp_regions.txt $params.sepsis_snps | awk '{ print \$16; }' | sort | uniq > snps_lists/sepsis_up_snps.txt
        """
}

process PREPARE_SNP_LIST_IN_LD {

    publishDir "$params.output_dir/snp_lists/", mode: "copy"

    label "simple_bash"

    input:
        path(snp_list)

    output:
        path("${snp_list.getSimpleName()}_ld.txt"), emit: snps_list
    
    script:

        """
        # Merge LD information from lead and conditional analyses
        awk 'FNR > 1 { print \$0; }' $params.lead_snps $params.conditional_snps | sort -k 1,1 | uniq > snp_ld_info.txt

        # Append SNPs that tag the lead SNPs with R^2 > 0.8
        join -j 1 -t \$'\\t' $snp_list snp_ld_info.txt | sed 's/NONE//g' | awk '
            BEGIN { OFS="\\t"; }
            {
                print NR, \$1;
                n_tags = split(\$8, tags, "|");
                for (i = 1; i <= n_tags; i++) {
                    print NR, tags[i];
                }
            }
        ' | sort -k 2,2 > snps_with_tags.txt

        # Merge with SNP information

        sort -k 2,2 $params.genotype_info > sorted_geno.bim

        echo -e "Locus\tSNP\tChr\tPosition" > ${snp_list.getSimpleName()}_ld.txt
        join -j 2 -t \$'\\t' sorted_geno.bim snps_with_tags.txt | awk -F "\\t" 'OFS="\\t" { print \$7, \$1, \$2, \$4; }' | sort -k 2,2n -k 3,3n >> ${snp_list.getSimpleName()}_ld.txt
        """
}

process PREPARE_PEAKS {

    publishDir "$params.output_dir/", mode: "copy"

    label "multi_cpu_bash"

    output:
        path("peaks/*.tsv"), emit: da_peak_sets
    
    script:

        """
        mkdir peaks/

        sed 's/"//g' $params.da_peaks | awk -F "," 'NR > 1 { print \$3; }' | sort | uniq > cell_types.txt

        while read cell_type; do

            echo -e "Chr\tStart\tEnd" >> peaks/\${cell_type}.tsv
        
            sed 's/"//g' $params.da_peaks | \\
                awk -F "," -v c=\$cell_type 'BEGIN { OFS="\\t"; } NR > 1 { if (\$3 == c) { print \$5, \$6, \$7; } }' | \\
                sort -k 1,1 -k 2,2n -k 3,3n >> peaks/\${cell_type}.tsv

        done <cell_types.txt
        """
}

process GO_SHIFTER {

    publishDir "$params.output_dir/", mode: "move"

    label "goshifter"

    input:
        path(snps_to_test)
        path(genome_annotation)

    output:
        path("${snps_to_test.getSimpleName()}_${genome_annotation.getSimpleName()}_p_value.tsv")
        path("${snps_to_test.getSimpleName()}_${genome_annotation.getSimpleName()}_overlap_scores.tsv")
    
    script:

        """
        go_shifter=\$(which lean_go_shifter.py)

        python \$go_shifter \\
            $snps_to_test \\
            $genome_annotation \\
            10000 \\
            ./ \\
            ${snps_to_test.getSimpleName()}_${genome_annotation.getSimpleName()} \\
            --threads $task.cpus
        """
}


workflow {

    PREPARE_SNP_LIST()

    PREPARE_SNP_LIST_IN_LD(PREPARE_SNP_LIST.out.snps_lists.flatten())

    PREPARE_PEAKS()

    snps_annotation_pairs = PREPARE_SNP_LIST_IN_LD.out.snps_list
        .flatten()
        .combine(PREPARE_PEAKS.out.da_peak_sets.flatten())
        .multiMap{ items ->
            snps_to_test: items[0]
            genome_annotation: items[1]
        }

    GO_SHIFTER(
        snps_annotation_pairs.snps_to_test,
        snps_annotation_pairs.genome_annotation
    )
}
