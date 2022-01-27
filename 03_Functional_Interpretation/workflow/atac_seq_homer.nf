nextflow.enable.dsl = 2

params.da_peak_set = "/nfs/users/nfs_n/nm18/gains_team282/epigenetics/calderon_et_al/analysis/atac_seq/da_peak_set.csv"
params.consensus_peak_set = "/nfs/users/nfs_n/nm18/gains_team282/epigenetics/calderon_et_al/analysis/atac_seq/consensus_peaks.bed"
params.genome_annotation = "/lustre/scratch118/humgen/resources/rna_seq_genomes/Homo_sapiens.GRCh38.99.gtf"
params.genome = "/lustre/scratch118/humgen/resources/rna_seq_genomes/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
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

    label "simple_bash"

    output:
        path("*.up.peaks.bed"),     emit: cell_specific_up_peaks
        path("*.down.peaks.bed"),   emit: cell_specific_down_peaks
    
    script:
    
        """
        awk -F ',' 'NR > 1 { gsub("\\"", "", \$3); print \$3; }' $params.da_peak_set | sort | uniq > cell_types.txt

        while read cell_type
        do

            sed 's/"//g' $params.da_peak_set | \\
                awk -F ',' -v c=\$cell_type 'NR > 1, OFS="\t" { if (\$3 == c && \$8 > 0) { print \$5, \$6, \$7, \$4, \$8, "."; } }' > \${cell_type}.up.peaks.bed
            
            sed 's/"//g' $params.da_peak_set | \\
                awk -F ',' -v c=\$cell_type 'NR > 1, OFS="\t" { if (\$3 == c && \$8 < 0) { print \$5, \$6, \$7, \$4, \$8, "."; } }' > \${cell_type}.down.peaks.bed
        
        done <cell_types.txt
        """
}

process PEAK_MOTIF_ENRICHMENT {

    errorStrategy "retry"
    maxRetries 3

    label "homer"

    publishDir "$params.output_dir/homer/${up_peak_set.getSimpleName()}/", mode: "move"

    input:
        path(up_peak_set)
        path(down_peak_set)

    output:
        path("*.known_motifs.tsv")
        path("*.all_motifs.motifs")
        path("*.peaks.tsv")

    script:

        """
        findMotifsGenome.pl \\
            $up_peak_set \\
            $params.genome \\
            ./homer_up_motifs/ \\
            -size given \\
            -p $task.cpus \\
            -preparsedDir ./preparsed/
        
        findMotifsGenome.pl \\
            $down_peak_set \\
            $params.genome \\
            ./homer_down_motifs/ \\
            -size given \\
            -p $task.cpus \\
            -preparsedDir ./preparsed/

        ls -1 homer_up_motifs/homerResults/ | \\
            grep "motif[0-9]*.motif" | \\
            sed 's/^/homer_up_motifs\\/homerResults\\//g' | \\
            xargs cat > homer_up_motifs/homerUpResults.motif

        annotatePeaks.pl \\
            $up_peak_set \\
            $params.genome -gtf $params.genome_annotation \\
            -m ./homer_up_motifs/homerUpResults.motif \\
            -noann -nogene > ${up_peak_set.getSimpleName()}.up.peaks.tsv
        
        ls -1 homer_down_motifs/homerResults/ | \\
            grep "motif[0-9]*.motif" | \\
            sed 's/^/homer_down_motifs\\/homerResults\\//g' | \\
            xargs cat > homer_down_motifs/homerDownResults.motif

        annotatePeaks.pl \\
            $down_peak_set \\
            $params.genome -gtf $params.genome_annotation \\
            -m ./homer_down_motifs/homerDownResults.motif \\
            -noann -nogene > ${down_peak_set.getSimpleName()}.down.peaks.tsv
        
        if [ -f "./homer_up_motifs/knownResults.txt" ];
        then
            mv ./homer_up_motifs/knownResults.txt ${up_peak_set.getSimpleName()}.up.known_motifs.tsv
        else
            touch ${up_peak_set.getSimpleName()}.up.known_motifs.tsv
        fi
        
        mv ./homer_up_motifs/homerMotifs.all.motifs ${up_peak_set.getSimpleName()}.up.all_motifs.motifs

        if [ -f "./homer_down_motifs/knownResults.txt" ];
        then
            mv ./homer_down_motifs/knownResults.txt ${down_peak_set.getSimpleName()}.down.known_motifs.tsv
        else
            touch ${down_peak_set.getSimpleName()}.down.known_motifs.tsv
        fi

        mv ./homer_down_motifs/homerMotifs.all.motifs ${down_peak_set.getSimpleName()}.down.all_motifs.motifs
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
