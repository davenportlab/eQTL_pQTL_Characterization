nextflow.enable.dsl = 2

params.da_peak_set = "/nfs/users/nfs_n/nm18/gains_team282/epigenetics/calderon_et_al/analysis/atac_seq/da_peak_set.csv"
params.consensus_peak_set = "/nfs/users/nfs_n/nm18/gains_team282/epigenetics/calderon_et_al/analysis/atac_seq/consensus_peaks.bed"
params.genome_annotation = "/lustre/scratch118/humgen/resources/rna_seq_genomes/Homo_sapiens.GRCh38.99.gtf"
params.genome = "/lustre/scratch118/humgen/resources/rna_seq_genomes/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
params.motifs = "/nfs/users/nfs_n/nm18/eQTL_pQTL_Characterization/03_Functional_Interpretation/data/JASPAR2022_CORE_vertebrates_non-redundant_pfms_meme.txt"
params.output_dir = "/nfs/users/nfs_n/nm18/gains_team282/epigenetics/calderon_et_al/analysis/atac_seq/"

//-----------------------------------------------------
// Annotate Peaks
//-----------------------------------------------------

process ANNOTATE_PEAKS {

    errorStrategy "retry"
    maxRetries 3

    label "homer"

    publishDir "$params.output_dir/homer/", mode: "move"
    
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

process DA_PEAKS_BY_CELL_TYPE {

    errorStrategy "retry"
    maxRetries 3

    label "bedtools"

    output:
        path("*.up.peaks.fa"),      emit: cell_specific_up_peaks
        path("*.down.peaks.fa"),    emit: cell_specific_down_peaks
    
    script:
    
        """
        cp $params.genome ./genome.fa
        samtools faidx genome.fa

        awk -F ',' 'NR > 1 { gsub("\\"", "", \$3); print \$3; }' $params.da_peak_set | sort | uniq > cell_types.txt

        while read cell_type
        do

            sed 's/"//g' $params.da_peak_set | \\
                awk -F ',' -v c=\$cell_type 'NR > 1, OFS="\t" { if (\$3 == c && \$8 > 0) { print \$5, \$6, \$7, \$4, \$8, "."; } }' > \${cell_type}.up.peaks.bed
            
            bedtools getfasta -fi genome.fa -bed \${cell_type}.up.peaks.bed > \${cell_type}.up.peaks.fa

            sed 's/"//g' $params.da_peak_set | \\
                awk -F ',' -v c=\$cell_type 'NR > 1, OFS="\t" { if (\$3 == c && \$8 < 0) { print \$5, \$6, \$7, \$4, \$8, "."; } }' > \${cell_type}.down.peaks.bed

            bedtools getfasta -fi genome.fa -bed \${cell_type}.down.peaks.bed > \${cell_type}.down.peaks.fa
        
        done <cell_types.txt
        """
}

process PEAK_MOTIF_ENRICHMENT {

    errorStrategy "retry"
    maxRetries 3

    label "sea"

    publishDir "$params.output_dir/homer/${up_peak_set.getSimpleName()}/", mode: "move"

    input:
        path(up_peak_set)
        path(down_peak_set)
        
    output:
        path("*.enrichment.tsv")
        path("*.motifs.tsv")
    
    script:

        """
        if [ -s $up_peak_set ] && [ \$(wc -l < $up_peak_set) -gt 2 ];
        then

            sea --seed 30023 --qvalue 0.05 --noseqs --p $up_peak_set --m $params.motifs --o up_motifs/
            
            mv up_motifs/sea.tsv ${up_peak_set.getSimpleName()}.up.enrichment.tsv

            fimo --o up_motifs_search/ $params.motifs $up_peak_set

            mv up_motifs_search/fimo.tsv ${up_peak_set.getSimpleName()}.up.motifs.tsv

        else

            touch ${up_peak_set.getSimpleName()}.up.enrichment.tsv
            touch ${up_peak_set.getSimpleName()}.up.motifs.tsv
        
        fi

        if [ -s $down_peak_set ] && [ \$(wc -l < $down_peak_set) -gt 2 ];
        then

            sea --seed 30023 --qvalue 0.05 --noseqs --p $down_peak_set --m $params.motifs --o down_motifs/

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

    DA_PEAKS_BY_CELL_TYPE()

    PEAK_MOTIF_ENRICHMENT(
        DA_PEAKS_BY_CELL_TYPE.out.cell_specific_up_peaks.flatten(),
        DA_PEAKS_BY_CELL_TYPE.out.cell_specific_down_peaks.flatten()
    )
}
