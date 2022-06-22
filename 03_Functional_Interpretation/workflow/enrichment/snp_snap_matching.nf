nextflow.enable.dsl = 2

params.geno_bim = "/nfs/users/nfs_n/nm18/gains_team282/Genotyping/All_genotyping_merged_filtered_b38_refiltered_rsID.bim"
params.immune_peak_sets = "/nfs/users/nfs_n/nm18/gains_team282/epigenetics/accessibility/analysis/atac_seq/immune/peak_sets/"
params.neutrophil_peak_sets = "/nfs/users/nfs_n/nm18/gains_team282/epigenetics/accessibility/analysis/atac_seq/neutrophil/peak_sets/"
params.encode_cres = "/nfs/users/nfs_n/nm18/gains_team282/nikhil/data/ENCODE/"
params.chromhmm = "/nfs/users/nfs_n/nm18/gains_team282/nikhil/data/ROADMAP/"
params.lead_snps = "/nfs/users/nfs_n/nm18/gains_team282/nikhil/colocalization/cis_eqtl/fine_mapping/LD/lead_snps.80r2.tags.tsv"
params.conditional_snps = "/nfs/users/nfs_n/nm18/gains_team282/nikhil/colocalization/cis_eqtl/fine_mapping/LD/conditional_snps.80r2.tags.tsv"
params.sepsis_snps = "/lustre/scratch119/humgen/projects/gains_team282/eqtl/gtex/sepsis_specificity_categories.txt"
params.snp_map = "/nfs/users/nfs_n/nm18/gains_team282/nikhil/data/SNPsnap/input_snps_identifier_mapping.txt"
params.matched_snps = "/nfs/users/nfs_n/nm18/gains_team282/nikhil/data/SNPsnap/matched_snps.txt"
params.hg19_to_hg38_chain = "/nfs/users/nfs_n/nm18/eQTL_pQTL_Characterization/03_Functional_Interpretation/data/hg19ToHg38.over.chain"
params.output_dir = "/nfs/users/nfs_n/nm18/gains_team282/epigenetics/enrichment/snp_snap/"


process GENOMIC_ANNOTATIONS {

    label "simple_bash"

    output:
        path("*.annotation.bed"),                   emit: annotations
        path("*.annotation.bed"),                   emit: annotations_for_matched
        path("lead_snps.bed"),                      emit: lead_snps
        path("conditional_snps.bed"),               emit: conditional_snps
        path("sepsis_enhanced_snps.bed"),           emit: sepsis_snps
        path("lead_matched_snps.txt"),              emit: lead_matched_snps
        path("conditional_matched_snps.txt"),       emit: conditional_matched_snps
        path("sepsis_enhanced_matched_snps.txt"),   emit: sepsis_matched_snps

    script:

        """
        ### Sort Genotype Information
        sort -k 2,2 $params.geno_bim > geno.bim

        ### Immune Atlas

        ls -1 $params.immune_peak_sets | grep \\.bed | sed 's/\\.peaks\\.bed//g' > immune_sets_names.txt

        while read immune_set_name
        do
            sort -k1,1 -k2,2n $params.immune_peak_sets/\${immune_set_name}.peaks.bed > ./\${immune_set_name}.annotation.bed
        done <immune_sets_names.txt

        ### Neutrophil Atlas

        ls -1 $params.neutrophil_peak_sets | grep \\.bed | sed 's/\\.peaks\\.bed//g' > neutrophil_sets_names.txt

        while read neutrophil_set_name
        do
            sort -k1,1 -k2,2n $params.neutrophil_peak_sets/\${neutrophil_set_name}.peaks.bed > ./\${neutrophil_set_name}.annotation.bed
        done <neutrophil_sets_names.txt

        ### ENCODE cCREs

        ls -1 $params.encode_cres | grep "GRCh38-cCREs\\\\..*\\\\.bed" | sed 's/\\.bed//g' > encode_cres_types.txt

        while read encode_cre_type
        do
            sed 's/^chr//g' $params.encode_cres/\${encode_cre_type}.bed | sort -k1,1 -k2,2n > ./\${encode_cre_type}.annotation.bed
        done <encode_cres_types.txt

        ### CHROM HMM

        ls -1 $params.chromhmm | grep '\\-E.*bed' | sed 's/\\.bed//g' > chromhmm_tracks.txt

        while read chromhmm_track
        do
            sed 's/^chr//g' $params.chromhmm/\${chromhmm_track}.bed | sort -k1,1 -k2,2n > ./\${chromhmm_track}.annotation.bed
        done <chromhmm_tracks.txt

        ### Make SNP Maps for cis-eQTL, conditional cis-eQTL, sepsis-enhanced SNPs, and sepsis-enhanced SNPs stronger in GAinS

        awk 'NR > 1 { print \$2 "\\t" \$1; }' $params.snp_map | sort -k1,1 > snp_map.txt
        
        # cis-eQTL
        awk 'NR > 1 { print \$1; }' $params.lead_snps | sort | join - snp_map.txt -t \$'\\t' | awk '{ print \$2; }' | sort > lead_snps.txt
        awk 'BEGIN { OFS="\\t"; } NR > 1 { print \$2, \$3, \$3 + 1; }' $params.lead_snps > lead_snps.bed

        # conditional cis-eQTL
        awk 'NR > 1 { print \$1; }' $params.conditional_snps | sort | join - snp_map.txt -t \$'\\t' | awk '{ print \$2; }' | sort > conditional_snps.txt
        awk 'BEGIN { OFS="\\t"; } NR > 1 { print \$2, \$3, \$3 + 1; }' $params.conditional_snps > conditional_snps.bed

        # sepsis-enhanced SNPs
        awk -F "\\t" '{ if (\$8 == "Bigger effect in GAinS") { print \$11; } }' $params.sepsis_snps | sort | join - snp_map.txt -t \$'\\t' | awk '{ print \$2; }' | sort > sepsis_enhanced_snps.txt
        awk -F "\\t" '{ if (\$8 == "Bigger effect in GAinS") { print \$11; } }' $params.sepsis_snps | sort | join -1 1 -2 2 - geno.bim -t \$'\\t' | awk 'OFS="\\t" { print \$2, \$4, \$4 + 1; }' > sepsis_enhanced_snps.bed

        ### Sort Matched SNPs and separate into cis-eQTL and conditional cis-eQTL

        awk 'NR > 1 { print \$0; }' $params.matched_snps | sort -k1,1 | join - lead_snps.txt > lead_matched_snps.txt

        awk 'NR > 1 { print \$0; }' $params.matched_snps | sort -k1,1 | join - conditional_snps.txt > conditional_matched_snps.txt

        awk 'NR > 1 { print \$0; }' $params.matched_snps | sort -k1,1 | join - sepsis_enhanced_snps.txt > sepsis_enhanced_matched_snps.txt
        """
}

process OVERLAP_OBSERVED {

    label "simple_bash"

    input:
        val(sets)
        path("annotations/*")
        path(lead_snps)
        path(conditional_snps)
        path(sepsis_snps)

    output:
        path("cis_eqtl_sets_results.tsv"),              emit: cis_eqtl_overlaps
        path("conditional_cis_eqtl_sets_results.tsv"),  emit: conditional_cis_eqtl_overlaps
        path("sepsis_enhanced_sets_results.tsv"),       emit: sepsis_overlaps

    script:

        """
        for i in ${sets.join(" ")}
        do

            ### cis-eQTL

            n_snps=\$(wc -l < $lead_snps)
            shuf -r -n \$n_snps $lead_snps | sort -k1,1 -k2,2n > set_hg38_sorted.bed

            for annotation in annotations/*
            do

                annotation_name=\$(basename \$annotation | sed 's/\\.annotation\\.bed//g')

                overlap=\$(bedtools intersect -a set_hg38_sorted.bed -b \$annotation -wa | wc -l)

                echo -e "\$i\\t\$annotation_name\\t\$overlap\\t\$n_snps" >> cis_eqtl_sets_results.tsv
            
            done

            ### Conditional cis-eQTL

            n_snps=\$(wc -l < $conditional_snps)
            shuf -r -n \$n_snps $conditional_snps | sort -k1,1 -k2,2n > set_hg38_sorted.bed

            for annotation in annotations/*
            do

                annotation_name=\$(basename \$annotation | sed 's/\\.annotation\\.bed//g')

                overlap=\$(bedtools intersect -a set_hg38_sorted.bed -b \$annotation -wa | wc -l)

                echo -e "\$i\\t\$annotation_name\\t\$overlap\\t\$n_snps" >> conditional_cis_eqtl_sets_results.tsv
            
            done

            ### Sepsis Enhanced cis-eQTL

            n_snps=\$(wc -l < $sepsis_snps)
            shuf -r -n \$n_snps $sepsis_snps | sort -k1,1 -k2,2n > set_hg38_sorted.bed

            for annotation in annotations/*
            do

                annotation_name=\$(basename \$annotation | sed 's/\\.annotation\\.bed//g')

                overlap=\$(bedtools intersect -a set_hg38_sorted.bed -b \$annotation -wa | wc -l)

                echo -e "\$i\\t\$annotation_name\\t\$overlap\\t\$n_snps" >> sepsis_enhanced_sets_results.tsv
            
            done

        done
        """
}

process LIFT_OVER_AND_OVERLAP {

    label "simple_bash"

    input:
        val(sets)
        path("annotations/*")
        path(lead_matched_snps)
        path(conditional_matched_snps)
        path(sepsis_matched_snps)
    
    output:
        path("cis_eqtl_sets_results.tsv"),              emit: cis_eqtl_overlaps
        path("conditional_cis_eqtl_sets_results.tsv"),  emit: conditional_cis_eqtl_overlaps
        path("sepsis_enhanced_sets_results.tsv"),       emit: sepsis_overlaps
    
    script:

        def set_indices = sets.collect{ it + 1 }

        """
        for i in ${set_indices.join(" ")}
        do

            ### cis-eQTL

            awk -v i=\$i '{ print \$i; }' $lead_matched_snps | sed 's/:/\\t/g' | awk 'OFS="\\t" { print "chr" \$1, \$2, \$2 + 1; }' > set.bed

            liftOver set.bed $params.hg19_to_hg38_chain set_hg38.bed unlifted.bed
            sed 's/^chr//g' set_hg38.bed | sort -k1,1 -k2,2n > set_hg38_sorted.bed

            for annotation in annotations/*
            do

                annotation_name=\$(basename \$annotation | sed 's/\\.annotation\\.bed//g')

                overlap=\$(bedtools intersect -a set_hg38_sorted.bed -b \$annotation -wa | wc -l)
                n_snps=\$(wc -l set_hg38_sorted.bed | sed 's/ set_hg38_sorted\\.bed//g')

                echo -e "\$i\\t\$annotation_name\\t\$overlap\\t\$n_snps" >> cis_eqtl_sets_results.tsv
            
            done

            ### Conditional cis-eQTL

            awk -v i=\$i '{ print \$i; }' $conditional_matched_snps | sed 's/:/\\t/g' | awk 'OFS="\\t" { print "chr" \$1, \$2, \$2 + 1; }' > set.bed

            liftOver set.bed $params.hg19_to_hg38_chain set_hg38.bed unlifted.bed
            sed 's/^chr//g' set_hg38.bed | sort -k1,1 -k2,2n > set_hg38_sorted.bed

            for annotation in annotations/*
            do

                annotation_name=\$(basename \$annotation | sed 's/\\.annotation\\.bed//g')

                overlap=\$(bedtools intersect -a set_hg38_sorted.bed -b \$annotation -wa | wc -l)
                n_snps=\$(wc -l set_hg38_sorted.bed | sed 's/ set_hg38_sorted\\.bed//g')

                echo -e "\$i\\t\$annotation_name\\t\$overlap\\t\$n_snps" >> conditional_cis_eqtl_sets_results.tsv
            
            done

            ### Sepsis Enhanced cis-eQTL

            awk -v i=\$i '{ print \$i; }' $sepsis_matched_snps | sed 's/:/\\t/g' | awk 'OFS="\\t" { print "chr" \$1, \$2, \$2 + 1; }' > set.bed

            liftOver set.bed $params.hg19_to_hg38_chain set_hg38.bed unlifted.bed
            sed 's/^chr//g' set_hg38.bed | sort -k1,1 -k2,2n > set_hg38_sorted.bed

            for annotation in annotations/*
            do

                annotation_name=\$(basename \$annotation | sed 's/\\.annotation\\.bed//g')

                overlap=\$(bedtools intersect -a set_hg38_sorted.bed -b \$annotation -wa | wc -l)
                n_snps=\$(wc -l set_hg38_sorted.bed | sed 's/ set_hg38_sorted\\.bed//g')

                echo -e "\$i\\t\$annotation_name\\t\$overlap\\t\$n_snps" >> sepsis_enhanced_sets_results.tsv
            
            done

        done
        """
}

process COLLATE_RESULTS {

    publishDir "$params.output_dir/", mode: "move"

    label "simple_bash"

    input:
        path("observed_cis_eqtl_results/*.tsv")
        path("observed_conditional_cis_eqtl_results/*.tsv")
        path("observed_sepsis_results/*.tsv")
        path("cis_eqtl_results/*.tsv")
        path("conditional_cis_eqtl_results/*.tsv")
        path("sepsis_results/*.tsv")
    
    output:
        path("*.tsv")

    script:

        """
        ### OBSERVED

        echo -e "Set\\tAnnotation\\tOverlap\\tN_SNPS" > cis_eqtl_observed_snp_overlaps.tsv
        cat observed_cis_eqtl_results/*.tsv | sort -k1,1n -k2,2 | awk 'OFS="\\t" { print \$1, \$2, \$3, \$4; }' >> cis_eqtl_observed_snp_overlaps.tsv

        echo -e "Set\\tAnnotation\\tOverlap\\tN_SNPS" > conditional_cis_eqtl_observed_snp_overlaps.tsv
        cat observed_conditional_cis_eqtl_results/*.tsv | sort -k1,1n -k2,2 | awk 'OFS="\\t" { print \$1, \$2, \$3, \$4; }' >> conditional_cis_eqtl_observed_snp_overlaps.tsv

        echo -e "Set\\tAnnotation\\tOverlap\\tN_SNPS" > sepsis_enhanced_observed_snp_overlaps.tsv
        cat observed_sepsis_results/*.tsv | sort -k1,1n -k2,2 | awk 'OFS="\\t" { print \$1, \$2, \$3, \$4; }' >> sepsis_enhanced_observed_snp_overlaps.tsv

        ### MATCHED

        echo -e "Set\\tAnnotation\\tOverlap\\tN_SNPS" > cis_eqtl_matched_snp_overlaps.tsv
        cat cis_eqtl_results/*.tsv | sort -k1,1n -k2,2 | awk 'OFS="\\t" { print \$1 - 1, \$2, \$3, \$4; }' >> cis_eqtl_matched_snp_overlaps.tsv

        echo -e "Set\\tAnnotation\\tOverlap\\tN_SNPS" > conditional_cis_eqtl_matched_snp_overlaps.tsv
        cat conditional_cis_eqtl_results/*.tsv | sort -k1,1n -k2,2 | awk 'OFS="\\t" { print \$1 - 1, \$2, \$3, \$4; }' >> conditional_cis_eqtl_matched_snp_overlaps.tsv

        echo -e "Set\\tAnnotation\\tOverlap\\tN_SNPS" > sepsis_enhanced_matched_snp_overlaps.tsv
        cat sepsis_results/*.tsv | sort -k1,1n -k2,2 | awk 'OFS="\\t" { print \$1 - 1, \$2, \$3, \$4; }' >> sepsis_enhanced_matched_snp_overlaps.tsv
        """
}

workflow {

    GENOMIC_ANNOTATIONS()

    OVERLAP_OBSERVED(
        Channel.from(1..10000).buffer(size: 100),
        GENOMIC_ANNOTATIONS.out.annotations.collect(),
        GENOMIC_ANNOTATIONS.out.lead_snps,
        GENOMIC_ANNOTATIONS.out.conditional_snps,
        GENOMIC_ANNOTATIONS.out.sepsis_snps
    )

    LIFT_OVER_AND_OVERLAP(
        Channel.from(1..10000).buffer(size: 100),
        GENOMIC_ANNOTATIONS.out.annotations_for_matched.collect(),
        GENOMIC_ANNOTATIONS.out.lead_matched_snps,
        GENOMIC_ANNOTATIONS.out.conditional_matched_snps,
        GENOMIC_ANNOTATIONS.out.sepsis_matched_snps
    )

    COLLATE_RESULTS(
        OVERLAP_OBSERVED.out.cis_eqtl_overlaps.collect(),
        OVERLAP_OBSERVED.out.conditional_cis_eqtl_overlaps.collect(),
        OVERLAP_OBSERVED.out.sepsis_overlaps.collect(),
        LIFT_OVER_AND_OVERLAP.out.cis_eqtl_overlaps.collect(),
        LIFT_OVER_AND_OVERLAP.out.conditional_cis_eqtl_overlaps.collect(),
        LIFT_OVER_AND_OVERLAP.out.sepsis_overlaps.collect()
    )
}
