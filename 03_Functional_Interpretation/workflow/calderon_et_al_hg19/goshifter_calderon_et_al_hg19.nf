nextflow.enable.dsl = 2

params.lead_snps_hg19 = "/nfs/users/nfs_n/nm18/gains_team282/epigenetics/calderon_et_al_hg19/lead_and_tag_snps_hg19.tsv"
params.conditional_snps_hg19 = "/nfs/users/nfs_n/nm18/gains_team282/epigenetics/calderon_et_al_hg19/conditional_and_tag_snps_hg19.tsv"
params.tss_regions_hg19 = "/nfs/users/nfs_n/nm18/gains_team282/epigenetics/calderon_et_al_hg19/tss_regions_hg19.tsv"
params.lead_snps = "/nfs/users/nfs_n/nm18/gains_team282/nikhil/colocalization/cis_eqtl/fine_mapping/LD/lead_snps.80r2.tags.tsv"
params.conditional_snps = "/nfs/users/nfs_n/nm18/gains_team282/nikhil/colocalization/cis_eqtl/fine_mapping/LD/conditional_snps.80r2.tags.tsv"
params.conditional_results = "/nfs/users/nfs_n/nm18/gains_team282/eqtl/cisresults/conditionalanalysis/conditional_eQTL_results_final.txt"
params.sepsis_snps = "/lustre/scratch119/humgen/projects/gains_team282/eqtl/gtex/sepsis_specific_gains_lead_mashr.txt"
params.mashr_results = "/lustre/scratch119/humgen/projects/gains_team282/eqtl/gtex/gains_gtex_mashr_results.txt"
params.da_peaks = "/nfs/users/nfs_n/nm18/eQTL_pQTL_Characterization/03_Functional_Interpretation/data/41588_2019_505_MOESM6_ESM"
params.goshifter_dir = "/nfs/users/nfs_n/nm18/gains_team282/nikhil/bin/goshifter/"
params.output_dir = "/nfs/users/nfs_n/nm18/gains_team282/epigenetics/calderon_et_al_hg19/"


process POSITIVE_CONTROL_EQTL_STUDIES {

    label "R"

    output:
        path("snps_lists/*.txt"), emit: snps_lists
    
    script:

        """
        mkdir snps_lists/

        Rscript $workflow.projectDir/calderon_et_al_hg19_positive_controls.R
        """
}

process PREPARE_SNP_LIST {

    label "simple_bash"

    output:
        path("snps_lists/*.txt"), emit: snps_lists
    
    script:

        """
        mkdir snps_lists/

        # Sort genotype information for merging in next steps
        awk 'NR > 1 { print \$0; }' $params.lead_snps_hg19 | awk 'OFS="\t" { print \$2, \$3, \$4; }' | sort -k 1,1 | uniq > sorted_lead_positions.tsv
        awk 'NR > 1 { print \$0; }' $params.conditional_snps_hg19 | awk 'OFS="\t" { print \$2, \$3, \$4; }' | sort -k 1,1 | uniq > sorted_conditional_positions.tsv

        ### Lead SNPs

        # Extract lead SNPs identified from initial cis-eQTL pass
        awk -F "\t" 'NR > 1 { print \$1; }' $params.lead_snps > lead_snp_IDs.txt

        # Merge with SNP information
        sort lead_snp_IDs.txt | uniq > sorted_lead_snps.txt
        join -1 1 -2 1 -t \$'\t' sorted_lead_positions.tsv sorted_lead_snps.txt | awk -F "\t" 'OFS="\t" { print \$1, "chr" \$2, \$3; }' > snps_lists/lead_snps.txt

        ### Lead SNPs from Conditional Analysis

        # Extract lead SNPs identified from conditional cis-eQTL pass
        awk -F "\t" 'NR > 1 { print \$1; }' $params.conditional_snps > conditional_snp_IDs.txt

        # Merge with SNP information
        sort conditional_snp_IDs.txt | uniq > sorted_conditional_snps.txt
        join -1 1 -2 1 -t \$'\t' sorted_conditional_positions.tsv sorted_conditional_snps.txt | awk -F "\t" 'OFS="\t" { print \$1, "chr" \$2, \$3; }' > snps_lists/conditional_snps.txt

        ### Lead SNPs from Conditional Analysis w/o Primary Effect

        # Extract lead SNPs identified from conditional cis-eQTL pass that are not the primary effect
        sed 's/"//g' $params.conditional_results | awk 'NR > 1 { if (\$7 > 1) { print \$2; } }' | sort | uniq > conditional_snp_IDs.txt

        # Merge with SNP information
        sort conditional_snp_IDs.txt | uniq > sorted_conditional_snps.txt
        join -1 1 -2 1 -t \$'\t' sorted_conditional_positions.tsv sorted_conditional_snps.txt | awk -F "\t" 'OFS="\t" { print \$1, "chr" \$2, \$3; }' > snps_lists/conditional_secondary_snps.txt

        ### Lead SNPs w/o Sepsis-Specific eSNPs

        # Extract lead SNPs identified from initial cis-eQTL pass
        awk -F "\t" 'NR > 1 { print \$1; }' $params.lead_snps > lead_snp_IDs.txt

        # Extract SNPs from mashR results
        awk 'NR > 1 { print \$16; }' $params.sepsis_snps | sort | uniq > sepsis_snps.txt

        # Filter lead SNPs to non-sepsis-specific eSNPs
        grep -vwFf sepsis_snps.txt lead_snp_IDs.txt | sort | uniq > non_sepsis_snps.txt

        # Merge with SNP information
        join -1 1 -2 1 -t \$'\t' sorted_lead_positions.tsv non_sepsis_snps.txt | awk -F "\t" 'OFS="\t" { print \$1, "chr" \$2, \$3; }' > snps_lists/non_sepsis_snps.txt

        ### Sepsis-Specific eSNPs

        # Extract SNPs from mashR results
        awk 'NR > 1 { print \$16; }' $params.sepsis_snps | sort | uniq > sepsis_snps.txt

        # Merge with SNP information
        join -1 1 -2 1 -t \$'\t' sorted_lead_positions.tsv sepsis_snps.txt | awk -F "\t" 'OFS="\t" { print \$1, "chr" \$2, \$3; }' > snps_lists/sepsis_snps.txt

        ### Sepsis-Specific eSNPs stronger in GAinS

        # Get SNP regions that are stronger in sepsis from the mashR results
        awk 'NR > 1 { if (\$5 == "TRUE" && (\$2)^2 > (\$3)^2) { print \$1; } }' $params.mashr_results | sort > sepsis_up_snp_regions.txt

        # Filter sepsis snps based on mashR results
        grep -wFf sepsis_up_snp_regions.txt $params.sepsis_snps | awk '{ print \$16; }' | sort | uniq > sepsis_up_snps.txt

        # Merge with SNP information
        join -1 1 -2 1 -t \$'\t' sorted_lead_positions.tsv sepsis_up_snps.txt | awk -F "\t" 'OFS="\t" { print \$1, "chr" \$2, \$3; }' > snps_lists/sepsis_up_snps.txt
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
        cp $snp_list snps_lists/$snp_list

        # Sort genotype information for merging in next steps
        awk 'NR > 1 { print \$0; }' $params.lead_snps_hg19 $params.conditional_snps_hg19 | awk 'OFS="\t" { print \$2, \$3, \$4; }' | sort -k 1,1 | uniq > sorted_positions.tsv

        # Merge LD information from lead and conditional analyses
        awk 'NR > 1 { print \$0; }' $params.lead_snps $params.conditional_snps | sort -k 1,1 | uniq > snp_ld_info.txt

        # Extract only the SNPs from the SNP list
        awk '{ print \$1; }' $snp_list | sort | uniq > snps.txt

        # Append SNPs that tag the lead SNPs with R^2 > 0.8
        join -1 1 -2 1 -t \$'\t' snps.txt snp_ld_info.txt | awk -F "\t" '{ if (\$4 > 0) { print \$8; } }' | sed 's/|/\\n/g' | sort | uniq > snps_with_tags.txt

        # Merge with SNP information
        join -1 1 -2 1 -t \$'\t' sorted_positions.tsv snps_with_tags.txt | awk -F "\t" 'OFS="\t" { print \$1, "chr" \$2, \$3; }' > snps_lists/${snp_list.getSimpleName()}_ld.txt
        """
}

process DA_PEAK_ANNOTATION_SETS {

    label "simple_bash"

    output:
        path("annotations/*.bed.gz")
    
    script:

        """
        mkdir annotations/

        awk 'NR > 1 { print \$8; }' $params.da_peaks | sort | uniq > contrasts.txt

        while read contrast
        do
            awk -v c="\$contrast" 'NR > 1 { if (\$8 == c) { gsub(/_/, "\t", \$7); print \$7; } }' $params.da_peaks > annotations/\${contrast}.bed
            gzip annotations/\${contrast}.bed
        done <contrasts.txt
        """
}

process GO_SHIFTER {

    publishDir "$params.output_dir/", mode: "move"

    label "goshifter"

    input:
        path(snps_to_test)

    output:
        path("goshifter/*.txt")
        path("goshifter/*.log")

    script:

        """
        mkdir goshifter/

        echo -e "SNP\tChrom\tBP" > snps.txt
        cat $snps_to_test >> snps.txt



        python $params.goshifter_dir/goshifter.py 
        """
}


workflow {

    POSITIVE_CONTROL_EQTL_STUDIES()

    PREPARE_SNP_LIST()

    PREPARE_SNP_LIST_IN_LD(PREPARE_SNP_LIST.out.snps_lists.flatten())

    all_snp_lists = PREPARE_SNP_LIST_IN_LD.out.snps_lists
        .concat(POSITIVE_CONTROL_EQTL_STUDIES.out.snps_lists)
        .flatten()

    GO_SHIFTER(all_snp_lists)
}
