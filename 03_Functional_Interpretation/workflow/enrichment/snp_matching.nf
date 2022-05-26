nextflow.enable.dsl = 2

params.output_dir = "/nfs/users/nfs_n/nm18/gains_team282/epigenetics/enrichment/snp_matching/"


process ENRICHMENT {

    publishDir "$params.output_dir/", mode: "move"

    label "multi_cpu_bash"

    output:
        path("fisher_results.csv")
        path("gains_tss_corrected_results.csv")
        path("shared_tss_corrected_results.csv")

    script:

        """
        Rscript $workflow.projectDir/snp_matching_generate_intervals.R

        sort -k1,1 -k2,2n cis_esnps.bed > cis_esnps.sorted.bed
        sort -k1,1 -k2,2n gains_esnps.bed > gains_esnps.sorted.bed
        sort -k1,1 -k2,2n shared_esnps.bed > shared_esnps.sorted.bed

        ## IMMUNE ATLAS

        ls -1 /nfs/users/nfs_n/nm18/gains_team282/epigenetics/accessibility/analysis/atac_seq/immune/peak_sets/ | grep "peaks.bed" > peak_set_beds.txt

        while read peak_set_bed
        do

            peak_type=\$(echo \$peak_set_bed | sed 's/\\.peaks\\.bed//g')

            bedtools intersect -a cis_esnps.sorted.bed -b /nfs/users/nfs_n/nm18/gains_team282/epigenetics/accessibility/analysis/atac_seq/immune/peak_sets/\${peak_set_bed} | awk '{ print \$4; }' | sort | uniq > \${peak_type}_snp_overlaps.txt
            awk '{ print \$4; }' cis_esnps.sorted.bed | sort | uniq | grep -vFf \${peak_type}_snp_overlaps.txt | awk '{ print \$0 "\\t0"; }' > cis_\${peak_type}_overlaps.txt
            awk '{ print \$0 "\\t1"; }' \${peak_type}_snp_overlaps.txt >> cis_\${peak_type}_overlaps.txt
            echo -e "Gene_SNP_Pair\\t\${peak_type}" > cis_\${peak_type}.overlaps.txt
            sort -k 1,1 cis_\${peak_type}_overlaps.txt >> cis_\${peak_type}.overlaps.txt

            bedtools intersect -a gains_esnps.sorted.bed -b /nfs/users/nfs_n/nm18/gains_team282/epigenetics/accessibility/analysis/atac_seq/immune/peak_sets/\${peak_set_bed} | awk '{ print \$4; }' | sort | uniq > \${peak_type}_snp_overlaps.txt
            awk '{ print \$4; }' gains_esnps.sorted.bed | sort | uniq | grep -vFf \${peak_type}_snp_overlaps.txt | awk '{ print \$0 "\\t0"; }' > gains_\${peak_type}_overlaps.txt
            awk '{ print \$0 "\\t1"; }' \${peak_type}_snp_overlaps.txt >> gains_\${peak_type}_overlaps.txt
            echo -e "Gene_SNP_Pair\\t\${peak_type}" > gains_\${peak_type}.overlaps.txt
            sort -k 1,1 gains_\${peak_type}_overlaps.txt >> gains_\${peak_type}.overlaps.txt

            bedtools intersect -a shared_esnps.sorted.bed -b /nfs/users/nfs_n/nm18/gains_team282/epigenetics/accessibility/analysis/atac_seq/immune/peak_sets/\${peak_set_bed} | awk '{ print \$4; }' | sort | uniq > \${peak_type}_snp_overlaps.txt
            awk '{ print \$4; }' shared_esnps.sorted.bed | sort | uniq | grep -vFf \${peak_type}_snp_overlaps.txt | awk '{ print \$0 "\\t0"; }' > shared_\${peak_type}_overlaps.txt
            awk '{ print \$0 "\\t1"; }' \${peak_type}_snp_overlaps.txt >> shared_\${peak_type}_overlaps.txt
            echo -e "Gene_SNP_Pair\\t\${peak_type}" > shared_\${peak_type}.overlaps.txt
            sort -k 1,1 shared_\${peak_type}_overlaps.txt >> shared_\${peak_type}.overlaps.txt

            rm \${peak_type}_snp_overlaps.txt cis_\${peak_type}_overlaps.txt gains_\${peak_type}_overlaps.txt shared_\${peak_type}_overlaps.txt

        done <peak_set_beds.txt

        ## NEUTROPHIL ATLAS

        ls -1 /nfs/users/nfs_n/nm18/gains_team282/epigenetics/accessibility/analysis/atac_seq/neutrophil/peak_sets/ | grep "peaks.bed" > peak_set_beds.txt

        while read peak_set_bed
        do

            peak_type=\$(echo \$peak_set_bed | sed 's/\\.peaks\\.bed//g')

            bedtools intersect -a cis_esnps.sorted.bed -b /nfs/users/nfs_n/nm18/gains_team282/epigenetics/accessibility/analysis/atac_seq/neutrophil/peak_sets/\${peak_set_bed} | awk '{ print \$4; }' | sort | uniq > \${peak_type}_snp_overlaps.txt
            awk '{ print \$4; }' cis_esnps.sorted.bed | sort | uniq | grep -vFf \${peak_type}_snp_overlaps.txt | awk '{ print \$0 "\\t0"; }' > cis_\${peak_type}_overlaps.txt
            awk '{ print \$0 "\\t1"; }' \${peak_type}_snp_overlaps.txt >> cis_\${peak_type}_overlaps.txt
            echo -e "Gene_SNP_Pair\\t\${peak_type}" > cis_\${peak_type}.overlaps.txt
            sort -k 1,1 cis_\${peak_type}_overlaps.txt >> cis_\${peak_type}.overlaps.txt

            bedtools intersect -a gains_esnps.sorted.bed -b /nfs/users/nfs_n/nm18/gains_team282/epigenetics/accessibility/analysis/atac_seq/neutrophil/peak_sets/\${peak_set_bed} | awk '{ print \$4; }' | sort | uniq > \${peak_type}_snp_overlaps.txt
            awk '{ print \$4; }' gains_esnps.sorted.bed | sort | uniq | grep -vFf \${peak_type}_snp_overlaps.txt | awk '{ print \$0 "\\t0"; }' > gains_\${peak_type}_overlaps.txt
            awk '{ print \$0 "\\t1"; }' \${peak_type}_snp_overlaps.txt >> gains_\${peak_type}_overlaps.txt
            echo -e "Gene_SNP_Pair\\t\${peak_type}" > gains_\${peak_type}.overlaps.txt
            sort -k 1,1 gains_\${peak_type}_overlaps.txt >> gains_\${peak_type}.overlaps.txt

            bedtools intersect -a shared_esnps.sorted.bed -b /nfs/users/nfs_n/nm18/gains_team282/epigenetics/accessibility/analysis/atac_seq/neutrophil/peak_sets/\${peak_set_bed} | awk '{ print \$4; }' | sort | uniq > \${peak_type}_snp_overlaps.txt
            awk '{ print \$4; }' shared_esnps.sorted.bed | sort | uniq | grep -vFf \${peak_type}_snp_overlaps.txt | awk '{ print \$0 "\\t0"; }' > shared_\${peak_type}_overlaps.txt
            awk '{ print \$0 "\\t1"; }' \${peak_type}_snp_overlaps.txt >> shared_\${peak_type}_overlaps.txt
            echo -e "Gene_SNP_Pair\\t\${peak_type}" > shared_\${peak_type}.overlaps.txt
            sort -k 1,1 shared_\${peak_type}_overlaps.txt >> shared_\${peak_type}.overlaps.txt

            rm \${peak_type}_snp_overlaps.txt cis_\${peak_type}_overlaps.txt gains_\${peak_type}_overlaps.txt shared_\${peak_type}_overlaps.txt

        done <peak_set_beds.txt

        ## ENCODE cCREs

        ls -1 /nfs/users/nfs_n/nm18/gains_team282/nikhil/data/ENCODE | grep "GRCh38-cCREs\\\\..*\\\\.bed" > peak_set_beds.txt

        while read peak_set_bed
        do

            peak_type=\$(echo \$peak_set_bed | sed 's/GRCh38-cCREs\\.//g' | sed 's/\\.bed//g')

            bedtools intersect -a cis_esnps.sorted.bed -b /nfs/users/nfs_n/nm18/gains_team282/nikhil/data/ENCODE/\${peak_set_bed} | awk '{ print \$4; }' | sort | uniq > \${peak_type}_snp_overlaps.txt
            awk '{ print \$4; }' cis_esnps.sorted.bed | sort | uniq | grep -vFf \${peak_type}_snp_overlaps.txt | awk '{ print \$0 "\\t0"; }' > cis_\${peak_type}_overlaps.txt
            awk '{ print \$0 "\\t1"; }' \${peak_type}_snp_overlaps.txt >> cis_\${peak_type}_overlaps.txt
            echo -e "Gene_SNP_Pair\\t\${peak_type}" > cis_\${peak_type}.overlaps.txt
            sort -k 1,1 cis_\${peak_type}_overlaps.txt >> cis_\${peak_type}.overlaps.txt

            bedtools intersect -a gains_esnps.sorted.bed -b /nfs/users/nfs_n/nm18/gains_team282/nikhil/data/ENCODE/\${peak_set_bed} | awk '{ print \$4; }' | sort | uniq > \${peak_type}_snp_overlaps.txt
            awk '{ print \$4; }' gains_esnps.sorted.bed | sort | uniq | grep -vFf \${peak_type}_snp_overlaps.txt | awk '{ print \$0 "\\t0"; }' > gains_\${peak_type}_overlaps.txt
            awk '{ print \$0 "\\t1"; }' \${peak_type}_snp_overlaps.txt >> gains_\${peak_type}_overlaps.txt
            echo -e "Gene_SNP_Pair\\t\${peak_type}" > gains_\${peak_type}.overlaps.txt
            sort -k 1,1 gains_\${peak_type}_overlaps.txt >> gains_\${peak_type}.overlaps.txt

            bedtools intersect -a shared_esnps.sorted.bed -b /nfs/users/nfs_n/nm18/gains_team282/nikhil/data/ENCODE/\${peak_set_bed} | awk '{ print \$4; }' | sort | uniq > \${peak_type}_snp_overlaps.txt
            awk '{ print \$4; }' shared_esnps.sorted.bed | sort | uniq | grep -vFf \${peak_type}_snp_overlaps.txt | awk '{ print \$0 "\\t0"; }' > shared_\${peak_type}_overlaps.txt
            awk '{ print \$0 "\\t1"; }' \${peak_type}_snp_overlaps.txt >> shared_\${peak_type}_overlaps.txt
            echo -e "Gene_SNP_Pair\\t\${peak_type}" > shared_\${peak_type}.overlaps.txt
            sort -k 1,1 shared_\${peak_type}_overlaps.txt >> shared_\${peak_type}.overlaps.txt

            rm \${peak_type}_snp_overlaps.txt cis_\${peak_type}_overlaps.txt gains_\${peak_type}_overlaps.txt shared_\${peak_type}_overlaps.txt

        done <peak_set_beds.txt

        paste cis_*.overlaps.txt | awk '{ printf "%s\t", \$1; for (i = 2; i <= NF; i += 2) { printf "%s", \$i; if (i != NF) { printf "\\t"; } else { printf "\\n"; } } }' > cis_overlaps.tsv
        paste gains_*.overlaps.txt | awk '{ printf "%s\t", \$1; for (i = 2; i <= NF; i += 2) { printf "%s", \$i; if (i != NF) { printf "\\t"; } else { printf "\\n"; } } }' > gains_overlaps.tsv
        paste shared_*.overlaps.txt | awk '{ printf "%s\t", \$1; for (i = 2; i <= NF; i += 2) { printf "%s", \$i; if (i != NF) { printf "\\t"; } else { printf "\\n"; } } }' > shared_overlaps.tsv
        
        rm cis_*.overlaps.txt
        rm gains_*.overlaps.txt
        rm shared_*.overlaps.txt

        Rscript $workflow.projectDir/snp_matching_tests.R
        """
}

workflow {

    ENRICHMENT()
}
