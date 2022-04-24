nextflow.enable.dsl = 2

params.da_peak_ia_set = "/nfs/users/nfs_n/nm18/gains_team282/epigenetics/accessibility/analysis/atac_seq/da_peak_ia_set.csv"
params.da_peak_rm_set = "/nfs/users/nfs_n/nm18/gains_team282/epigenetics/accessibility/analysis/atac_seq/da_peak_rm_set.csv"
params.consensus_peak_set = "/nfs/users/nfs_n/nm18/gains_team282/epigenetics/accessibility/analysis/atac_seq/consensus_peaks.bed"
params.genome_annotation = "/lustre/scratch118/humgen/resources/rna_seq_genomes/Homo_sapiens.GRCh38.99.gtf"
params.genome = "/lustre/scratch118/humgen/resources/rna_seq_genomes/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
params.motifs = "/nfs/users/nfs_n/nm18/eQTL_pQTL_Characterization/03_Functional_Interpretation/data/JASPAR2022_CORE_vertebrates_non-redundant_pfms_meme.txt"
params.output_homer_dir = "/nfs/users/nfs_n/nm18/gains_team282/epigenetics/accessibility/analysis/atac_seq/homer/"
params.output_ia_dir = "/nfs/users/nfs_n/nm18/gains_team282/epigenetics/accessibility/analysis/atac_seq/immune_atlas_motifs/"
params.output_na_dir = "/nfs/users/nfs_n/nm18/gains_team282/epigenetics/accessibility/analysis/atac_seq/neutrophil_atlas_motifs/"

//-----------------------------------------------------
// Annotate Peaks
//-----------------------------------------------------

process ANNOTATE_PEAKS {

    errorStrategy "retry"
    maxRetries 3

    label "homer"

    publishDir "$params.output_homer_dir/", mode: "move"
    
    output:
        path("homer_peak_annotation.txt")

    script:

        """
        awk -F '\t' 'OFS="\t" { print \$1, \$2, \$3, \$1 ":" \$2 "-" \$3, \$4, "."; }' \\
            $params.consensus_peak_set > consensus_peaks_homer.bed

        annotatePeaks.pl \\
            consensus_peaks_homer.bed \\
            $params.genome -gtf $params.genome_annotation \\
            -organism human > homer_peak_annotation.txt
        """
}

//-----------------------------------------------------
// Motif Enrichment
//-----------------------------------------------------

process DA_PEAK_SEQUENCES {

    errorStrategy "retry"
    maxRetries 3

    label "bedtools"

    output:
        path("control.sequences.fa"),   emit: control_sequences
        path("*.ia.up.peaks.fa"),       emit: ia_up_peaks
        path("*.ia.down.peaks.fa"),     emit: ia_down_peaks
        path("*.rm.up.peaks.fa"),       emit: rm_up_peaks
        path("*.rm.down.peaks.fa"),     emit: rm_down_peaks
    
    script:
    
        """
        cp $params.genome ./genome.fa
        samtools faidx genome.fa

        # Control sequences, which are non-DA peaks
        sed 's/"//g' $params.da_peak_ia_set | \\
            awk -F ',' 'NR > 1, OFS="\t" { print \$5, \$6, \$7, \$4, "0", "."; }' > da.peaks.ia.bed
        
        sed 's/"//g' $params.da_peak_rm_set | \\
            awk -F ',' 'NR > 1, OFS="\t" { print \$5, \$6, \$7, \$4, "0", "."; }' > da.peaks.rm.bed

        cat da.peaks.ia.bed da.peaks.rm.bed | sort | uniq > da.peaks.bed
        
        bedtools subtract -a $params.consensus_peak_set -b da.peaks.bed -A > control.sequences.bed

        bedtools getfasta -fi genome.fa -bed control.sequences.bed > control.sequences.fa

        # Determine contrasts
        #   For immune atlas, this is the cell type
        #   For neutrophil atlas, this is the treatment
        awk -F ',' 'NR > 1 { gsub("\\"", "", \$3); print \$3; }' $params.da_peak_ia_set | sort | uniq > cell_types.txt
        awk -F ',' 'NR > 1 { gsub("\\"", "", \$3); print \$3; }' $params.da_peak_rm_set | sort | uniq > treatment_types.txt
        
        while read cell_type
        do

            sed 's/"//g' $params.da_peak_ia_set | \\
                awk -F ',' -v c=\$cell_type 'NR > 1, OFS="\t" { if (\$3 == c && \$8 > 0) { print \$5, \$6, \$7, \$4, \$8, "."; } }' > \${cell_type}.up.peaks.bed
            
            bedtools getfasta -fi genome.fa -bed \${cell_type}.up.peaks.bed > \${cell_type}.ia.up.peaks.fa

            sed 's/"//g' $params.da_peak_ia_set | \\
                awk -F ',' -v c=\$cell_type 'NR > 1, OFS="\t" { if (\$3 == c && \$8 < 0) { print \$5, \$6, \$7, \$4, \$8, "."; } }' > \${cell_type}.down.peaks.bed

            bedtools getfasta -fi genome.fa -bed \${cell_type}.down.peaks.bed > \${cell_type}.ia.down.peaks.fa
        
        done <cell_types.txt

        while read treatment
        do

            sed 's/"//g' $params.da_peak_rm_set | \\
                awk -F ',' -v c=\$treatment 'NR > 1, OFS="\t" { if (\$3 == c && \$8 > 0) { print \$5, \$6, \$7, \$4, \$8, "."; } }' > \${treatment}.up.peaks.bed
            
            bedtools getfasta -fi genome.fa -bed \${treatment}.up.peaks.bed > \${treatment}.rm.up.peaks.fa

            sed 's/"//g' $params.da_peak_rm_set | \\
                awk -F ',' -v c=\$treatment 'NR > 1, OFS="\t" { if (\$3 == c && \$8 < 0) { print \$5, \$6, \$7, \$4, \$8, "."; } }' > \${treatment}.down.peaks.bed

            bedtools getfasta -fi genome.fa -bed \${treatment}.down.peaks.bed > \${treatment}.rm.down.peaks.fa
        
        done <treatment_types.txt
        """
}

process PEAK_MOTIF_ENRICHMENT_IA_UP {

    errorStrategy "retry"
    maxRetries 3

    label "sea"

    publishDir "$params.output_ia_dir/${up_peak_set.getSimpleName()}/", mode: "move"

    input:
        path(control_sequences)
        path(up_peak_set)
        
    output:
        path("*.enrichment.tsv")
        path("*.motifs.tsv")
    
    script:

        """
        if [ -s $up_peak_set ] && [ \$(wc -l < $up_peak_set) -gt 2 ];
        then

            sea --seed 30023 --qvalue 0.05 --noseqs --p $up_peak_set --m $params.motifs --n $control_sequences --o up_motifs/
            
            mv up_motifs/sea.tsv ${up_peak_set.getSimpleName()}.up.enrichment.tsv

            fimo --o up_motifs_search/ $params.motifs $up_peak_set

            mv up_motifs_search/fimo.tsv ${up_peak_set.getSimpleName()}.up.motifs.tsv

        else

            touch ${up_peak_set.getSimpleName()}.up.enrichment.tsv
            touch ${up_peak_set.getSimpleName()}.up.motifs.tsv
        
        fi
        """
}

process PEAK_MOTIF_ENRICHMENT_IA_DOWN {

    errorStrategy "retry"
    maxRetries 3

    label "sea"

    publishDir "$params.output_ia_dir/${down_peak_set.getSimpleName()}/", mode: "move"

    input:
        path(control_sequences)
        path(down_peak_set)
        
    output:
        path("*.enrichment.tsv")
        path("*.motifs.tsv")
    
    script:

        """
        if [ -s $down_peak_set ] && [ \$(wc -l < $down_peak_set) -gt 2 ];
        then

            sea --seed 30023 --qvalue 0.05 --noseqs --p $down_peak_set --m $params.motifs --n $control_sequences --o down_motifs/

            mv down_motifs/sea.tsv ${down_peak_set.getSimpleName()}.down.enrichment.tsv

            fimo --o down_motifs_search/ $params.motifs $down_peak_set

            mv down_motifs_search/fimo.tsv ${down_peak_set.getSimpleName()}.down.motifs.tsv
        
        else

            touch ${down_peak_set.getSimpleName()}.down.enrichment.tsv
            touch ${down_peak_set.getSimpleName()}.down.motifs.tsv
        
        fi
        """
}

process PEAK_MOTIF_ENRICHMENT_RM_UP {

    errorStrategy "retry"
    maxRetries 3

    label "sea"

    publishDir "$params.output_na_dir/${up_peak_set.getSimpleName()}/", mode: "move"

    input:
        path(control_sequences)
        path(up_peak_set)
        
    output:
        path("*.enrichment.tsv")
        path("*.motifs.tsv")
    
    script:

        """
        if [ -s $up_peak_set ] && [ \$(wc -l < $up_peak_set) -gt 2 ];
        then

            sea --seed 30023 --qvalue 0.05 --noseqs --p $up_peak_set --m $params.motifs --n $control_sequences --o up_motifs/
            
            mv up_motifs/sea.tsv ${up_peak_set.getSimpleName()}.up.enrichment.tsv

            fimo --o up_motifs_search/ $params.motifs $up_peak_set

            mv up_motifs_search/fimo.tsv ${up_peak_set.getSimpleName()}.up.motifs.tsv

        else

            touch ${up_peak_set.getSimpleName()}.up.enrichment.tsv
            touch ${up_peak_set.getSimpleName()}.up.motifs.tsv
        
        fi
        """
}

process PEAK_MOTIF_ENRICHMENT_RM_DOWN {

    errorStrategy "retry"
    maxRetries 3

    label "sea"

    publishDir "$params.output_na_dir/${down_peak_set.getSimpleName()}/", mode: "move"

    input:
        path(control_sequences)
        path(down_peak_set)
        
    output:
        path("*.enrichment.tsv")
        path("*.motifs.tsv")
    
    script:

        """
        if [ -s $down_peak_set ] && [ \$(wc -l < $down_peak_set) -gt 2 ];
        then

            sea --seed 30023 --qvalue 0.05 --noseqs --p $down_peak_set --m $params.motifs --n $control_sequences --o down_motifs/

            mv down_motifs/sea.tsv ${down_peak_set.getSimpleName()}.down.enrichment.tsv

            fimo --o down_motifs_search/ $params.motifs $down_peak_set

            mv down_motifs_search/fimo.tsv ${down_peak_set.getSimpleName()}.down.motifs.tsv
        
        else

            touch ${down_peak_set.getSimpleName()}.down.enrichment.tsv
            touch ${down_peak_set.getSimpleName()}.down.motifs.tsv
        
        fi
        """
}

//-----------------------------------------------------
// Workflow
//-----------------------------------------------------

workflow {

    ANNOTATE_PEAKS()

    DA_PEAK_SEQUENCES()

    // Immune Atlas
    PEAK_MOTIF_ENRICHMENT_IA_UP(
        DA_PEAK_SEQUENCES.out.control_sequences,
        DA_PEAK_SEQUENCES.out.ia_up_peaks.flatten()
    )

    PEAK_MOTIF_ENRICHMENT_IA_DOWN(
        DA_PEAK_SEQUENCES.out.control_sequences,
        DA_PEAK_SEQUENCES.out.ia_down_peaks.flatten()
    )

    // Neutrophil Atlas
    PEAK_MOTIF_ENRICHMENT_RM_UP(
        DA_PEAK_SEQUENCES.out.control_sequences,
        DA_PEAK_SEQUENCES.out.rm_up_peaks.flatten()
    )

    PEAK_MOTIF_ENRICHMENT_RM_DOWN(
        DA_PEAK_SEQUENCES.out.control_sequences,
        DA_PEAK_SEQUENCES.out.rm_down_peaks.flatten()
    )
}
