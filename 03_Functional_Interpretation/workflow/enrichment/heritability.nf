nextflow.enable.dsl = 2

params.immune_peak_sets = "/nfs/users/nfs_n/nm18/gains_team282/epigenetics/accessibility/analysis/atac_seq/immune/peak_sets/"
params.neutrophil_peak_sets = "/nfs/users/nfs_n/nm18/gains_team282/epigenetics/accessibility/analysis/atac_seq/neutrophil/peak_sets/"
params.encode_cres = "/nfs/users/nfs_n/nm18/gains_team282/nikhil/data/ENCODE/"
params.chromhmm = "/nfs/users/nfs_n/nm18/gains_team282/nikhil/data/ROADMAP/"
params.genotypes_prefix = "/nfs/users/nfs_n/nm18/gains_team282/Genotyping/All_genotyping_merged_filtered_b38_refiltered_rsID"
params.mapping_patients = "/nfs/users/nfs_n/nm18/gains_team282/nikhil/expression/eigengene_sva/mapping_patients.txt"
params.output_dir = "/nfs/users/nfs_n/nm18/gains_team282/epigenetics/enrichment/heritability/"


process GENOMIC_ANNOTATIONS {

    label "simple_bash"

    output:
        path("*.annotation.bed"), emit: annotations

    script:
    
        """
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
        """
}

process GRM_SETUP {

    label "simple_bash"

    publishDir "$params.output_dir/", mode: "copy"

    output:
        path("autosomes.*"),    emit: autosome_snps
        path("full_grm.grm.*")

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
        """
}


process GRM_CALCULATION {

    publishDir "$params.output_dir/", mode: "move"

    input:
        path(annotation_bed)
        path("autosomes/*")

    output:
        path("${annotation_bed.getSimpleName()}/")

    script:

        """
        mkdir ${annotation_bed.getSimpleName()}/

        awk 'OFS="\\t" { print \$1, \$4 - 1, \$4, \$2; }' autosomes/autosomes.bim > snps.bed

        bedtools intersect -a snps.bed -b $annotation_bed -wa | awk '{ print \$4; }' > in_annotation.snps
        bedtools intersect -v -a snps.bed -b $annotation_bed -wa | awk '{ print \$4; }' > not_in_annotation.snps

        gcta64 \\
            --bfile autosomes/autosomes \\
            --extract in_annotation.snps \\
            --make-grm \\
            --out ${annotation_bed.getSimpleName()}/annotation_snps_grm \\
            --thread-num 16

        gcta64 \\
            --bfile autosomes/autosomes \\
            --extract not_in_annotation.snps \\
            --make-grm \\
            --out ${annotation_bed.getSimpleName()}/other_snps_grm \\
            --thread-num 16
        """
}

workflow {

    GENOMIC_ANNOTATIONS()

    GRM_SETUP()

    GRM_CALCULATION(
        GENOMIC_ANNOTATIONS.out.annotations.flatten(),
        GRM_SETUP.out.autosome_snps
    )
}
