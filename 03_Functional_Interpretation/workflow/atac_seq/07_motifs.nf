nextflow.enable.dsl = 2

params.atlas = "immune"
params.da_peak_set = "/nfs/users/nfs_n/nm18/gains_team282/epigenetics/accessibility/analysis/atac_seq/$params.atlas/da_peak_set.csv"
params.consensus_peak_set = "/nfs/users/nfs_n/nm18/gains_team282/epigenetics/accessibility/analysis/atac_seq/$params.atlas/consensus_peaks.bed"
params.cell_peak_sets_dir = "/nfs/users/nfs_n/nm18/gains_team282/epigenetics/accessibility/analysis/atac_seq/$params.atlas/cell_type_peak_sets/"
params.group_peak_sets_dir = "/nfs/users/nfs_n/nm18/gains_team282/epigenetics/accessibility/analysis/atac_seq/$params.atlas/peak_sets/"
params.genome_annotation = "/lustre/scratch118/humgen/resources/rna_seq_genomes/Homo_sapiens.GRCh38.99.gtf"
params.genome = "/lustre/scratch118/humgen/resources/rna_seq_genomes/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
params.chr_lengths = "/nfs/users/nfs_n/nm18/gains_team282/epigenetics/star_genome_index/chrNameLength.txt"
params.motifs = "/nfs/users/nfs_n/nm18/eQTL_pQTL_Characterization/03_Functional_Interpretation/data/JASPAR2022_CORE_vertebrates_non-redundant_pfms_meme.txt"
params.output_dir = "/nfs/users/nfs_n/nm18/gains_team282/epigenetics/accessibility/analysis/atac_seq/$params.atlas/motifs/"

//-----------------------------------------------------
// Annotate Peaks
//-----------------------------------------------------

process ANNOTATE_CONSENSUS_PEAKS {

    errorStrategy "retry"
    maxRetries 3

    label "homer"

    publishDir "$params.output_dir/homer/", mode: "move"
    
    output:
        path("consensus_peak_annotation.txt")

    script:

        """
        awk -F '\\t' 'OFS="\\t" { print \$1, \$2, \$3, \$1 ":" \$2 "-" \$3, \$4, "."; }' \\
            $params.consensus_peak_set > consensus_peaks_homer.bed

        annotatePeaks.pl \\
            consensus_peaks_homer.bed \\
            $params.genome -gtf $params.genome_annotation \\
            -organism human > consensus_peak_annotation.txt
        """
}

process ANNOTATE_CELL_PEAKS {

    errorStrategy "retry"
    maxRetries 3

    label "homer"

    publishDir "$params.output_dir/homer/cell_types/", mode: "move"
    
    output:
        path("*.peak_annotation.txt")

    script:

        """
        ls -1 $params.cell_peak_sets_dir | grep bed | sed 's/\\.peaks\\.bed//g' > cell_types.txt

        while read cell_type
        do

            awk -F '\\t' 'OFS="\\t" { print \$1, \$2, \$3, \$1 ":" \$2 "-" \$3, \$4, "."; }' \\
                $params.cell_peak_sets_dir/\${cell_type}.peaks.bed > peaks_homer.bed

            annotatePeaks.pl \\
                peaks_homer.bed \\
                $params.genome -gtf $params.genome_annotation \\
                -organism human > \${cell_type}.peak_annotation.txt
        
        done <cell_types.txt
        """
}

process ANNOTATE_GROUP_PEAKS {

    errorStrategy "retry"
    maxRetries 3

    label "homer"

    publishDir "$params.output_dir/homer/groups/", mode: "move"
    
    output:
        path("*.peak_annotation.txt")

    script:

        """
        ls -1 $params.group_peak_sets_dir | grep bed | sed 's/\\.peaks\\.bed//g' > groups.txt

        while read group
        do

            awk -F '\\t' 'OFS="\\t" { print \$1, \$2, \$3, \$1 ":" \$2 "-" \$3, \$4, "."; }' \\
                $params.group_peak_sets_dir/\${group}.peaks.bed > peaks_homer.bed

            annotatePeaks.pl \\
                peaks_homer.bed \\
                $params.genome -gtf $params.genome_annotation \\
                -organism human > \${group}.peak_annotation.txt
        
        done <groups.txt
        """
}

//-----------------------------------------------------
// Consensus Peak Set Motif Enrichment
//-----------------------------------------------------

process CONSENSUS_PEAK_SEQUENCES {

    errorStrategy "retry"
    maxRetries 3

    label "bedtools"

    output:
        path("control.sequences.fa"),   emit: control_sequences
        path("consensus.peaks.fa"),     emit: consensus_sequences

    script:

        """
        cp $params.genome ./genome.fa
        samtools faidx genome.fa

        bedtools flank \\
            -i $params.consensus_peak_set \\
            -g $params.chr_lengths \\
            -b 1.0 -pct | \\
            bedtools subtract -a - -b $params.consensus_peak_set -A > control.sequences.bed

        bedtools getfasta -fi genome.fa -bed control.sequences.bed > control.sequences.fa

        bedtools getfasta -fi genome.fa -bed $params.consensus_peak_set > consensus.peaks.fa
        """
}

process CONSENSUS_MOTIFS {

    errorStrategy "retry"
    maxRetries 3

    label "sea"

    publishDir "$params.output_dir/consensus/", mode: "move"

    input:
        path(control_sequences)
        path(peak_sequences)
        
    output:
        path("*.enrichment.tsv")
        path("*.motifs.tsv")
    
    script:

        """
        sea --seed 30023 --qvalue 0.05 --noseqs --p $peak_sequences --m $params.motifs --n $control_sequences --o motifs/
        
        mv motifs/sea.tsv consensus.enrichment.tsv

        fimo --o motifs_search/ $params.motifs $peak_sequences

        mv motifs_search/fimo.tsv consensus.motifs.tsv
        """
}

//-----------------------------------------------------
// Cell Type Peak Set Motif Enrichment
//-----------------------------------------------------

process CELL_PEAK_SEQUENCES {

    errorStrategy "retry"
    maxRetries 3

    label "bedtools"

    output:
        path("*.control.sequences.fa"), emit: control_sequences
        path("*.peaks.fa"),             emit: cell_sequences

    script:

        """
        cp $params.genome ./genome.fa
        samtools faidx genome.fa

        ls -1 $params.cell_peak_sets_dir | grep bed | sed 's/\\.peaks\\.bed//g' > cell_types.txt

        while read cell_type
        do

            bedtools flank \\
                -i $params.cell_peak_sets_dir/\${cell_type}.peaks.bed \\
                -g $params.chr_lengths \\
                -b 1.0 -pct | \\
                bedtools subtract -a - -b $params.cell_peak_sets_dir/\${cell_type}.peaks.bed -A > control.sequences.bed

            bedtools getfasta -fi genome.fa -bed control.sequences.bed > \${cell_type}.control.sequences.fa

            bedtools getfasta -fi genome.fa -bed $params.cell_peak_sets_dir/\${cell_type}.peaks.bed > \${cell_type}.peaks.fa

        done <cell_types.txt
        """
}

process CELL_MOTIFS {

    errorStrategy "retry"
    maxRetries 3

    label "sea"

    publishDir "$params.output_dir/cell_type_motifs/${peak_sequences.getSimpleName()}/", mode: "move"

    input:
        path(control_sequences)
        path(peak_sequences)
        
    output:
        path("*.enrichment.tsv")
        path("*.motifs.tsv")
    
    script:

        """
        if [ -s $peak_sequences ] && [ \$(wc -l < $peak_sequences) -gt 2 ];
        then

            sea --seed 30023 --qvalue 0.05 --noseqs --p $peak_sequences --m $params.motifs --n $control_sequences --o motifs/
            
            mv motifs/sea.tsv ${peak_sequences.getSimpleName()}.enrichment.tsv

            fimo --o motifs_search/ $params.motifs $peak_sequences

            mv motifs_search/fimo.tsv ${peak_sequences.getSimpleName()}.motifs.tsv
        
        else

            touch ${peak_sequences.getSimpleName()}.enrichment.tsv
            touch ${peak_sequences.getSimpleName()}.motifs.tsv
        
        fi
        """
}

//-----------------------------------------------------
// Group Peak Set Motif Enrichment
//-----------------------------------------------------

process GROUP_PEAK_SEQUENCES {

    errorStrategy "retry"
    maxRetries 3

    label "bedtools"

    output:
        path("*.control.sequences.fa"), emit: control_sequences
        path("*.peaks.fa"),             emit: group_sequences

    script:

    if (params.atlas == "immune")

        """
        cp $params.genome ./genome.fa
        samtools faidx genome.fa

        ls -1 $params.group_peak_sets_dir | grep bed | grep Control | sed 's/-Control\\.peaks\\.bed//g' | awk '{ print \$1 "\\t" \$1 "-Control"; }' > control_groups.txt
        ls -1 $params.group_peak_sets_dir | grep bed | grep Stimulated | sed 's/-Stimulated\\.peaks\\.bed//g' | awk '{ print \$1 "\\t" \$1 "-Stimulated"; }' > stimulated_groups.txt
        join control_groups.txt stimulated_groups.txt | awk '{ print \$1; }' > common_groups.txt

        while read group
        do

            bedtools intersect \\
                -a $params.group_peak_sets_dir/\${group}-Stimulated.peaks.bed \\
                -b $params.group_peak_sets_dir/\${group}-Control.peaks.bed \\
                -f 0.9 \\
                -wa -u > stimulated.overlaps.bed
            
            bedtools intersect \\
                -a $params.group_peak_sets_dir/\${group}-Control.peaks.bed \\
                -b $params.group_peak_sets_dir/\${group}-Stimulated.peaks.bed \\
                -f 0.9 \\
                -wa -u > control.overlaps.bed

            cat control.overlaps.bed stimulated.overlaps.bed | \\
                sort -k1,1 -k2,2n | \\
                bedtools merge -i - > control.sequences.bed

            bedtools subtract \\
                -a $params.group_peak_sets_dir/\${group}-Stimulated.peaks.bed \\
                -b $params.group_peak_sets_dir/\${group}-Control.peaks.bed \\
                -A > stimulated.sequences.bed

            bedtools getfasta -fi genome.fa -bed control.sequences.bed > \${group}.control.sequences.fa

            bedtools getfasta -fi genome.fa -bed stimulated.sequences.bed > \${group}.peaks.fa

        done <common_groups.txt
        """

    else if (params.atlas == "neutrophil")

        """
        cp $params.genome ./genome.fa
        samtools faidx genome.fa

        stimulated=( LTA LPS FLAG R848 BGP HMGB1 "SA-1" "SA-3" "SA-5" EC1h EC4h )
        controls=( Control Control Control Control Control Control WB WB WB noEC1h noEC4h )

        for i in {0..10}
        do

            bedtools intersect \\
                -a $params.group_peak_sets_dir/Neutrophils-\${stimulated[\$i]}.peaks.bed \\
                -b $params.group_peak_sets_dir/Neutrophils-\${controls[\$i]}.peaks.bed \\
                -f 0.9 \\
                -wa -u > stimulated.overlaps.bed
            
            bedtools intersect \\
                -a $params.group_peak_sets_dir/Neutrophils-\${controls[\$i]}.peaks.bed \\
                -b $params.group_peak_sets_dir/Neutrophils-\${stimulated[\$i]}.peaks.bed \\
                -f 0.9 \\
                -wa -u > control.overlaps.bed

            cat control.overlaps.bed stimulated.overlaps.bed | \\
                sort -k1,1 -k2,2n | \\
                bedtools merge -i - > control.sequences.bed

            bedtools subtract \\
                -a $params.group_peak_sets_dir/Neutrophils-\${stimulated[\$i]}.peaks.bed \\
                -b $params.group_peak_sets_dir/Neutrophils-\${controls[\$i]}.peaks.bed \\
                -A > stimulated.sequences.bed

            bedtools getfasta -fi genome.fa -bed control.sequences.bed > Neutrophils-\${stimulated[\$i]}.control.sequences.fa

            bedtools getfasta -fi genome.fa -bed stimulated.sequences.bed > Neutrophils-\${stimulated[\$i]}.peaks.fa

        done
        """
}

process GROUP_MOTIFS {

    errorStrategy "retry"
    maxRetries 3

    label "sea"

    publishDir "$params.output_dir/group_motifs/${peak_sequences.getSimpleName()}/", mode: "move"

    input:
        path(control_sequences)
        path(peak_sequences)
        
    output:
        path("*.enrichment.tsv")
        path("*.motifs.tsv")
    
    script:

        """
        if [ -s $peak_sequences ] && [ \$(wc -l < $peak_sequences) -gt 2 ];
        then

            sea --seed 30023 --qvalue 0.05 --noseqs --p $peak_sequences --m $params.motifs --n $control_sequences --o motifs/
            
            mv motifs/sea.tsv ${peak_sequences.getSimpleName()}.enrichment.tsv

            fimo --o motifs_search/ $params.motifs $peak_sequences

            mv motifs_search/fimo.tsv ${peak_sequences.getSimpleName()}.motifs.tsv
        
        else

            touch ${peak_sequences.getSimpleName()}.enrichment.tsv
            touch ${peak_sequences.getSimpleName()}.motifs.tsv
        
        fi
        """
}

//-----------------------------------------------------
// DA Peak Set Motif Enrichment
//-----------------------------------------------------

process DA_PEAK_SEQUENCES {

    errorStrategy "retry"
    maxRetries 3

    label "bedtools"

    output:
        path("control.sequences.fa"),   emit: control_sequences
        path("*.up.peaks.fa"),          emit: up_peaks
        path("*.down.peaks.fa"),        emit: down_peaks
    
    script:
    
        """
        cp $params.genome ./genome.fa
        samtools faidx genome.fa

        # Control sequences, which are non-DA peaks
        sed 's/"//g' $params.da_peak_set | \\
            awk -F ',' 'NR > 1, OFS="\\t" { print \$5, \$6, \$7, \$4, "0", "."; }' | \\
            sort | uniq > da.peaks.bed
        
        bedtools subtract -a $params.consensus_peak_set -b da.peaks.bed -A > control.sequences.bed

        bedtools getfasta -fi genome.fa -bed control.sequences.bed > control.sequences.fa

        # Determine contrasts
        #   For immune atlas, this is the cell type
        #   For neutrophil atlas, this is the treatment
        awk -F ',' 'NR > 1 { gsub("\\"", "", \$3); print \$3; }' $params.da_peak_set | sort | uniq > contrasts.txt
        
        while read contrast
        do

            sed 's/"//g' $params.da_peak_set | \\
                awk -F ',' -v c=\$contrast 'NR > 1, OFS="\\t" { if (\$3 == c && \$8 > 0) { print \$5, \$6, \$7, \$4, \$8, "."; } }' > \${contrast}.up.peaks.bed
            
            bedtools getfasta -fi genome.fa -bed \${contrast}.up.peaks.bed > \${contrast}.up.peaks.fa

            sed 's/"//g' $params.da_peak_set | \\
                awk -F ',' -v c=\$contrast 'NR > 1, OFS="\\t" { if (\$3 == c && \$8 < 0) { print \$5, \$6, \$7, \$4, \$8, "."; } }' > \${contrast}.down.peaks.bed

            bedtools getfasta -fi genome.fa -bed \${contrast}.down.peaks.bed > \${contrast}.down.peaks.fa
        
        done <contrasts.txt
        """
}

process DA_MOTIFS_UP {

    errorStrategy "retry"
    maxRetries 3

    label "sea"

    publishDir "$params.output_dir/da_peaks/${up_peak_set.getSimpleName()}/", mode: "move"

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

process DA_MOTIFS_DOWN {

    errorStrategy "retry"
    maxRetries 3

    label "sea"

    publishDir "$params.output_dir/da_peaks/${down_peak_set.getSimpleName()}/", mode: "move"

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

    ANNOTATE_CONSENSUS_PEAKS()

    ANNOTATE_CELL_PEAKS()

    ANNOTATE_GROUP_PEAKS()

    CONSENSUS_PEAK_SEQUENCES()

    CONSENSUS_MOTIFS(
        CONSENSUS_PEAK_SEQUENCES.out.control_sequences,
        CONSENSUS_PEAK_SEQUENCES.out.consensus_sequences
    )

    CELL_PEAK_SEQUENCES()

    CELL_MOTIFS(
        CELL_PEAK_SEQUENCES.out.control_sequences.flatten(),
        CELL_PEAK_SEQUENCES.out.cell_sequences.flatten()
    )

    GROUP_PEAK_SEQUENCES()

    GROUP_MOTIFS(
        GROUP_PEAK_SEQUENCES.out.control_sequences.flatten(),
        GROUP_PEAK_SEQUENCES.out.group_sequences.flatten()
    )

    DA_PEAK_SEQUENCES()

    DA_MOTIFS_UP(
        DA_PEAK_SEQUENCES.out.control_sequences,
        DA_PEAK_SEQUENCES.out.up_peaks.flatten()
    )

    DA_MOTIFS_DOWN(
        DA_PEAK_SEQUENCES.out.control_sequences,
        DA_PEAK_SEQUENCES.out.down_peaks.flatten()
    )
}
