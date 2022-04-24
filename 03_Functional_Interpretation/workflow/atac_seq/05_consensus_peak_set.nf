nextflow.enable.dsl = 2

params.atac_seq_dir = "/nfs/users/nfs_n/nm18/gains_team282/epigenetics/accessibility/merged/atac_seq/"
params.metadata = "/nfs/users/nfs_n/nm18/eQTL_pQTL_Characterization/03_Functional_Interpretation/metadata/reads_atac_seq.txt"
params.output_dir = "/nfs/users/nfs_n/nm18/gains_team282/epigenetics/accessibility/analysis/atac_seq/"


//-----------------------------------------------------
// Consensus Peak Set
//-----------------------------------------------------

process GROUP_PEAK_SET {

    errorStrategy "retry"
    maxRetries 3

    label "bedtools"

    publishDir "$params.output_dir/$atlas/peak_sets/", mode: "copy", pattern: "*.{bed,saf}"

    input:
        val(atlas)
        val(group)
        path("narrow_peaks/*.narrowPeak")

    output:
        tuple val(atlas), val(group), path("${group}.peaks.bed"), emit: group_peak_set
        tuple val(atlas), val(group), path("${group}.peaks.bed"), emit: group_peak_set_for_cell_type
        path("${group}.peaks.saf")

    script:

        // -d 0 - Maximum distance allowed is 0 / Overlap must be at least 1bp
        // -c 1,7,8,9,10 - Operate on columns 1 (Chromosome Name), 7 (Fold Enrichment), 8 (-log10(P-Value)), 9 (-log10(Q-Value)), and 10 (Point-Source for Peak)
        // -o count,collapse,collapse,collapse,collapse - Count peaks using column 1 and collapse columns 7,8,9,10 into a comma-separated list

        // Filtering via AWK
        //  1. Any peak with q-value >= 10E-4 is removed
        //  2. Any peak that is wider than 3 kb is removed
        //  3. Any peak only present in one sample is removed
        """
        cat narrow_peaks/*.narrowPeak | awk '{ if (\$9 > 4) { print \$0; } }' | sort -k1,1 -k2,2n > narrow_peaks.sorted.bed

        bedtools merge \\
            -d 0 \\
            -c 1,7,8,9,10 \\
            -o count,collapse,collapse,collapse,collapse \\
            -i narrow_peaks.sorted.bed | \\
            awk -F '\t' -v OFS='\t' '{
                if (\$3 - \$2 <= 3000) {
                    if (\$4 >= 2) {
                        print \$0;
                    }
                }
            }' | \\
            sort -k1,1 -k2,2n > ${group}.peaks.bed

        echo -e "GeneID\tChr\tStart\tEnd\tStrand" > ${group}.peaks.saf
        awk 'OFS="\t" { print \$1 ":" \$2 "-" \$3, \$1, \$2, \$3, "+"; }' ${group}.peaks.bed >> ${group}.peaks.saf
        """
}

process CELL_TYPE_PEAK_SET {

    errorStrategy "retry"
    maxRetries 3

    label "bedtools"

    publishDir "$params.output_dir/$atlas/cell_type_peak_sets/", mode: "copy", pattern: "*.{bed,saf}"

    input:
        val(atlas)
        val(cell_type)
        path("peak_sets/*")

    output:
        path("${cell_type}.peaks.bed")
        path("${cell_type}.peaks.saf")
    
    script:

        // Pre-Filtering
        //  1. A peak from a condition must overlap some other peak with 90% of the interval to qualify for merging
        //  2. A peak will also qualify for merging if it does not overlap a peak from any other condition

        // -d 0 - Maximum distance allowed is 0 / Overlap must be at least 1bp

        // Filtering via AWK
        //  1. Any peak that is wider than 3 kb is removed
        """
        mkdir peak_sets_pruned/

        num_peak_files=\$(ls -1 peak_sets/*.bed | wc -l)

        if [ \$num_peak_files -gt 1 ]
        then

            for peak_set in peak_sets/*.bed
            do

                other_peaks=\$(ls -1 peak_sets/*.bed | grep -v "\${peak_set}" | tr '\\n' ' ')
                bedtools intersect -a \$peak_set -b \$other_peaks -f 0.9 -wa -u > peak_sets_pruned/\$(basename \$peak_set)
                bedtools intersect -a \$peak_set -b \$other_peaks -wa -v >> peak_sets_pruned/\$(basename \$peak_set)

            done

        else

            cp peak_sets/*.bed peak_sets_pruned/

        fi

        cat peak_sets_pruned/*.bed | sort -k1,1 -k2,2n > peaks.sorted.bed

        bedtools merge \\
            -d 0 \\
            -i peaks.sorted.bed | \\
            awk 'OFS="\\t" { print \$1, \$2, \$3, \$1 ":" \$2 "-" \$3, "."; }' | \\
            awk -F '\t' -v OFS='\t' '{ if (\$3 - \$2 <= 3000) { print \$0; } }' | \\
            sort -k1,1 -k2,2n > ${cell_type}.peaks.bed

        echo -e "GeneID\tChr\tStart\tEnd\tStrand" > ${cell_type}.peaks.saf
        awk 'OFS="\t" { print \$4, \$1, \$2, \$3, "+"; }' ${cell_type}.peaks.bed >> ${cell_type}.peaks.saf
        """
}

process CONSENSUS_PEAK_SET {

    errorStrategy "retry"
    maxRetries 3

    label "bedtools"

    publishDir "$params.output_dir/$atlas/", mode: "copy", pattern: "*.{bed,saf}"

    input:
        val(atlas)
        path("peak_sets/*.bed")

    output:
        path("consensus_peaks.bed")
        path("consensus_peaks.saf")
    
    script:

        // Pre-Filtering
        //  1. A peak from a condition must overlap some other peak with 90% of the interval to qualify for merging
        //  2. A peak will also qualify for merging if it does not overlap a peak from any other condition

        // -d 0 - Maximum distance allowed is 0 / Overlap must be at least 1bp

        // Filtering via AWK
        //  1. Any peak that is wider than 3 kb is removed
        """
        mkdir peak_sets_pruned/

        for peak_set in peak_sets/*.bed
        do

            other_peaks=\$(ls -1 peak_sets/*.bed | grep -v "\${peak_set}" | tr '\\n' ' ')
            bedtools intersect -a \$peak_set -b \$other_peaks -f 0.9 -wa -u > peak_sets_pruned/\$(basename \$peak_set)
            bedtools intersect -a \$peak_set -b \$other_peaks -wa -v >> peak_sets_pruned/\$(basename \$peak_set)

        done

        cat peak_sets_pruned/*.bed | sort -k1,1 -k2,2n > peaks.sorted.bed

        bedtools merge \\
            -d 0 \\
            -i peaks.sorted.bed | \\
            awk 'OFS="\\t" { print \$1, \$2, \$3, \$1 ":" \$2 "-" \$3, "."; }' | \\
            awk -F '\t' -v OFS='\t' '{ if (\$3 - \$2 <= 3000) { print \$0; } }' | \\
            sort -k1,1 -k2,2n > consensus_peaks.bed

        echo -e "GeneID\tChr\tStart\tEnd\tStrand" > consensus_peaks.saf
        awk 'OFS="\t" { print \$4, \$1, \$2, \$3, "+"; }' consensus_peaks.bed >> consensus_peaks.saf
        """
}

//-----------------------------------------------------
// Workflow
//-----------------------------------------------------

workflow {

    //-----------------------------------------------------
    // Consensus Peak Set
    //-----------------------------------------------------

    immune_peak_sets = Channel
        .fromPath(params.metadata)
        .splitCsv(skip: 1)
        .filter{ row -> row[3] =~ /Corces/ || row[3] =~ /Calderon/ }
        .unique{ row -> row[1] }
        .map{ row -> [ [ row[5] + (row[6] == "TRUE" ? "-Stimulated" : "-Control"), "immune" ], file("$params.atac_seq_dir/${row[1]}/peaks/${row[1]}_peaks.narrowPeak") ] }

    neutrophil_peak_sets = Channel
        .fromPath(params.metadata)
        .splitCsv(skip: 1, sep: ",")
        .filter{ row -> row[3] =~ /Ram-Mohan/ }
        .unique{ row -> row[1] }
        .map{ row -> [ [ row[5] + "-" + row[7], "neutrophil"] , file("$params.atac_seq_dir/${row[1]}/peaks/${row[1]}_peaks.narrowPeak") ] }

     peak_sets = immune_peak_sets
        .concat(neutrophil_peak_sets)
        .groupTuple()
        .multiMap{ pair -> 
            group: pair[0][0]
            atlas: pair[0][1]
            peaks: pair[1]
        }

    GROUP_PEAK_SET(peak_sets.atlas, peak_sets.group, peak_sets.peaks)

    cell_type_peak_files = GROUP_PEAK_SET.out.group_peak_set_for_cell_type
        .map{ it -> [it[1].replaceAll(/-.*/, ""), it[0], it[2]] }
        .groupTuple()
        .multiMap{ pair ->
            cell_type: pair[0]
            atlas: pair[1][0]
            peaks: pair[2]
        }
    
    CELL_TYPE_PEAK_SET(
        cell_type_peak_files.atlas,
        cell_type_peak_files.cell_type,
        cell_type_peak_files.peaks
    )

    peak_files = GROUP_PEAK_SET.out.group_peak_set
        .groupTuple()
        .multiMap{ pair ->
            atlas: pair[0]
            peaks: pair[2]
        }

    CONSENSUS_PEAK_SET(peak_files.atlas, peak_files.peaks)
}
