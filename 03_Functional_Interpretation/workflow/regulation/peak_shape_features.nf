nextflow.enable.dsl = 2

params.atlas = "immune"
params.metadata = "/nfs/users/nfs_n/nm18/eQTL_pQTL_Characterization/03_Functional_Interpretation/metadata/reads_atac_seq.txt"
params.peak_dir = "/nfs/users/nfs_n/nm18/gains_team282/epigenetics/accessibility/analysis/atac_seq/immune/cell_type_peak_sets/"
params.consensus_peaks = "/nfs/users/nfs_n/nm18/gains_team282/epigenetics/accessibility/analysis/atac_seq/immune/consensus_peaks.bed"
params.output_dir = "/nfs/users/nfs_n/nm18/gains_team282/epigenetics/regulation/immune/"


//-----------------------------------------------------
// Calculate Shape Features
//-----------------------------------------------------

process SHAPE_FEATURES {

    errorStrategy "retry"
    maxRetries 3

    label "python"

    publishDir "$params.output_dir/shape_features/", mode: "move"

    input:
        path(peaks)
        val(cell_type)
    
    output:
        path("${cell_type}_shape_features.csv")
        path("${cell_type}_consensus_shape_features.csv")

    script:

        // I first filter the peaks to only autosomes and the X chromosome

        """
        grep -E "^[1-9]|^X" $peaks > cell_type_peaks.bed
        python3 $workflow.projectDir/peak_shape_features.py cell_type_peaks.bed $cell_type ${cell_type}_shape_features.csv

        grep -E "^[1-9]|^X" $params.consensus_peaks > consensus_peaks.bed
        python3 $workflow.projectDir/peak_shape_features.py consensus_peaks.bed $cell_type ${cell_type}_consensus_shape_features.csv
        """
}

//-----------------------------------------------------
// Workflow
//-----------------------------------------------------

workflow {

    samples_list = Channel
        .fromPath(params.metadata)
        .splitCsv(skip: 1)
    
    if (params.atlas == "immune") {
        samples_list = samples_list
            .filter{ row -> row[3] =~ /Corces/ || row[3] =~ /Calderon/ }
    }

    if (params.atlas == "neutrophil") {
        samples_list = samples_list
            .filter{ row -> row[3] =~ /Ram-Mohan/ }
    }

    samples_list = samples_list
        .unique{ row -> row[5] }
        .multiMap{ row -> 
            peaks: file("$params.peak_dir/${row[5]}.peaks.bed")
            cell_type: row[5]
        }

    SHAPE_FEATURES(samples_list.peaks, samples_list.cell_type)
}
