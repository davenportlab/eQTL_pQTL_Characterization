nextflow.enable.dsl = 2

params.genotype_info = "/nfs/users/nfs_n/nm18/gains_team282/Genotyping/All_genotyping_merged_filtered_b38_refiltered_rsID.bim"
params.gene_info = "/nfs/team282/data/gains_team282/gene_info_864_20412_hla.txt"
params.lead_snps = "/nfs/users/nfs_n/nm18/gains_team282/nikhil/colocalization/cis_eqtl/fine_mapping/LD/lead_snps.80r2.tags.tsv"
params.conditional_snps = "/nfs/users/nfs_n/nm18/gains_team282/nikhil/colocalization/cis_eqtl/fine_mapping/LD/conditional_snps.80r2.tags.tsv"
params.conditional_results = "/nfs/users/nfs_n/nm18/gains_team282/eqtl/cisresults/conditionalanalysis/conditional_eQTL_results_final.txt"
params.sepsis_snps = "/lustre/scratch119/humgen/projects/gains_team282/eqtl/gtex/sepsis_specific_gains_lead_mashr.txt"
params.mashr_results = "/lustre/scratch119/humgen/projects/gains_team282/eqtl/gtex/gains_gtex_mashr_results.txt"
params.peak_counts = "/nfs/users/nfs_n/nm18/gains_team282/epigenetics/accessibility/analysis/atac_seq/neutrophil/peak_counts.tsv"
params.output_dir = "/nfs/users/nfs_n/nm18/gains_team282/epigenetics/enrichment/"


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

        ### MERGE WITH SNP COORDINATES

        sort -k 2,2 $params.genotype_info > sorted_geno.bim
        
        for snp_list in snps_lists/*.txt
        do
            join -1 2 -2 1 -t \$'\\t' sorted_geno.bim \$snp_list | awk -F "\\t" 'OFS="\\t" { print \$1, "chr" \$2, \$4; }' > snps_lists/temp.txt
            mv snps_lists/temp.txt \$snp_list
        done
        """
}

process PREPARE_SNP_LIST_IN_LD {

    label "simple_bash"

    input:
        path(snp_list)

    output:
        path("snps_lists/*.txt"), emit: snps_lists
    
    script:

        """
        mkdir snps_lists/
        # cp $snp_list snps_lists/$snp_list

        # Merge LD information from lead and conditional analyses
        awk 'FNR > 1 { print \$0; }' $params.lead_snps $params.conditional_snps | sort -k 1,1 | uniq > snp_ld_info.txt

        # Extract only the SNPs from the SNP list
        awk '{ print \$1; }' $snp_list | sort | uniq > snps.txt

        # Append SNPs that tag the lead SNPs with R^2 > 0.8
        join -1 1 -2 1 -t \$'\\t' snps.txt snp_ld_info.txt | awk -F "\\t" '{ if (\$4 > 0) { print \$8; } }' | sed 's/|/\\n/g' | sort | uniq > snps_with_tags.txt

        # Merge with SNP information

        sort -k 2,2 $params.genotype_info > sorted_geno.bim

        join -1 2 -2 1 -t \$'\\t' sorted_geno.bim snps_with_tags.txt | awk -F "\\t" 'OFS="\\t" { print \$1, "chr" \$2, \$4; }' > snps_lists/${snp_list.getSimpleName()}_ld.txt
        """
}

process PREPARE_PEAKS {

    label "multi_cpu_bash"

    output:
        path("peaks/*.txt"), emit: merged_peak_counts
    
    script:

        """
        # Filter peak count data so that only peaks with at least one read in one sample are included
        
        awk '
            NR == 1 { print \$0; }
            NR > 1 {
                sum = 0;
                for (i = 6; i <= NF; i++) {
                    sum += \$i;
                }
                if (sum > 0) {
                    print \$0;
                }
            }' $params.peak_counts > peak_counts_filtered.txt

        # Filter peak count data so that only peaks in 1 Mb of an expressed gene's TSS are included
        
        cp $params.gene_info gene_info.txt
        sed 's/"//g' gene_info.txt | awk -F "\\t" '
            BEGIN { OFS="\\t"; } 
            NR > 1 { 
                if (\$2 != "NA") { 
                    start = 0;
                    end = 0;
                    if (\$6 == "+") { 
                        start = \$3 - 10^6;
                        end = \$3 + 10^6 - 1; 
                    } else { 
                        start = \$4 - 10^6;
                        end = \$4 + 10^6 - 1; 
                    }
                    if (start < 1) {
                        start = 1;
                    }
                    print "chr" \$2, start, end;
                } 
            }' > tss_regions.bed

        awk 'BEGIN { OFS="\\t"; } NR > 1 { print "chr" \$2, \$3, \$4; }' peak_counts_filtered.txt > filtered_peaks.bed

        bedtools intersect -a filtered_peaks.bed -b tss_regions.bed -wa -u | awk '{ gsub(/chr/, "", \$1); print \$1 ":" \$2 "-" \$3; }' > filtered_peaks_in_cis.txt

        head -n 1 peak_counts_filtered.txt > peak_counts_filtered_in_cis.txt
        grep -wFf filtered_peaks_in_cis.txt peak_counts_filtered.txt >> peak_counts_filtered_in_cis.txt

        # Calculate median peak width
        median_peak_width=\$(awk 'NR > 1 { print \$4 - \$3 + 1; }' peak_counts_filtered_in_cis.txt | datamash median 1)
        peak_window=\$(echo "\$median_peak_width / 2" | bc)

        # Merge peak counts for each cell type

        mkdir peaks/

        function process_cell_type {

            awk -v c="-\${1}" -v w="\${2}" '
                OFS="\\t",
                NR == 1 {
                    for (i = 1; i <= NF; i++) {
                        if (\$i ~ c) {
                            matching_columns[\$i] = i;
                        }
                    }
                }
                NR > 1 {
                    sum = 0;
                    for (key in matching_columns) {
                        sum += \$(matching_columns[key]);
                    }
                    peak_center = int((\$3 + \$4) / 2);
                    print "chr" \$2, peak_center - w, peak_center + w - 1, sum;
                }
                ' peak_counts_filtered_in_cis.txt > peaks/\${1}.txt
        }

        export -f process_cell_type

        parallel process_cell_type ::: LTA LPS FLAG R848 BGP HMGB1 Control WB SA-1 SA-3 SA-5 noEC1h noEC4h EC1h EC4h ::: \$peak_window
        """
}

process CHEERS_NORMALIZE {

    publishDir "$params.output_dir/cheers/", mode: "copy"

    label "cheers"

    input:
        path("peaks/*")
    
    output:
        path("normalized/*.txt"), emit: normalized_data

    script:

        """
        cheers_normalize=\$(which CHEERS_normalize.py)

        mkdir normalized/

        python \$cheers_normalize Neutrophil_Atlas normalized/ \\
            peaks/LTA.txt peaks/LPS.txt peaks/FLAG.txt peaks/R848.txt peaks/BGP.txt peaks/HMGB1.txt peaks/Control.txt

        python \$cheers_normalize Neutrophil_Atlas_WB normalized/ \\
            peaks/WB.txt peaks/SA-1.txt peaks/SA-3.txt peaks/SA-5.txt peaks/noEC1h.txt peaks/noEC4h.txt peaks/EC1h.txt peaks/EC4h.txt 
        """
}

process CHEERS {

    publishDir "$params.output_dir/", mode: "move"

    label "cheers"

    input:
        path(snps_to_test)
        path("normalized/*")

    output:
        path("cheers/*.txt")
        path("cheers/*.log")
    
    script:

        """
        cheers_enrichment=\$(which CHEERS_computeEnrichment.py)

        mkdir cheers/

        python \$cheers_enrichment \\
            ${snps_to_test.getSimpleName()} \\
            cheers/ \\
            normalized/Neutrophil_Atlas_counts_normToMax_quantileNorm_euclideanNorm.txt \\
            --snp_list $snps_to_test
        
        python \$cheers_enrichment \\
            ${snps_to_test.getSimpleName()}_WB \\
            cheers/ \\
            normalized/Neutrophil_Atlas_WB_counts_normToMax_quantileNorm_euclideanNorm.txt \\
            --snp_list $snps_to_test
        """
}


workflow {

    PREPARE_SNP_LIST()

    PREPARE_SNP_LIST_IN_LD(PREPARE_SNP_LIST.out.snps_lists.flatten())

    PREPARE_PEAKS()

    CHEERS_NORMALIZE(PREPARE_PEAKS.out.merged_peak_counts)

    CHEERS(
        PREPARE_SNP_LIST_IN_LD.out.snps_lists.flatten(),
        CHEERS_NORMALIZE.out.normalized_data
    )
}
