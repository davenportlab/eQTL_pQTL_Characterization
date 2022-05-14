nextflow.enable.dsl = 2

params.atlas = "immune"
params.metadata = "/nfs/users/nfs_n/nm18/eQTL_pQTL_Characterization/03_Functional_Interpretation/metadata/reads_atac_seq.txt"
params.peak_dir = "/nfs/users/nfs_n/nm18/gains_team282/epigenetics/accessibility/analysis/atac_seq/immune/cell_type_peak_sets/"
params.consensus_peaks = "/nfs/users/nfs_n/nm18/gains_team282/epigenetics/accessibility/analysis/atac_seq/immune/consensus_peaks.bed"
params.output_dir = "/nfs/users/nfs_n/nm18/gains_team282/epigenetics/regulation/immune/"


//-----------------------------------------------------
// Calculate Shape Features
//-----------------------------------------------------

process CONSENSUS_SHAPE_FEATURES {

    label "consensus_shape_features"

    publishDir "$params.output_dir/shape_features/", mode: "move"

    output:
        path("consensus_fiedler_vectors.csv")
        path("consensus_coverage_pcs.csv")
        path("consensus_coverage_pcs_var_explained.csv")
    
    script:

        def sample_regex = (params.atlas == "immune") ? "Corces|Calderon" : "Ram-Mohan"

        """
        awk -F ',' 'NR > 1 { print \$2; }' $params.metadata | grep -E "$sample_regex" | sort | uniq > samples.txt

        grep -E "^[1-9]|^X" $params.consensus_peaks > consensus_peaks.bed
        python3 $workflow.projectDir/peak_shape_features.py consensus_peaks.bed samples.txt consensus --threads 16
        """
}

process CELL_TYPE_SHAPE_FEATURES {

    label "cell_shape_features"

    publishDir "$params.output_dir/shape_features/", mode: "move"

    input:
        path(peaks)
        val(cell_type)
    
    output:
        path("${cell_type}_fiedler_vectors.csv")
        path("${cell_type}_coverage_pcs.csv")
        path("${cell_type}_coverage_pcs_var_explained.csv")

    script:

        // I first filter the peaks to only autosomes and the X chromosome

        """
        awk -F ',' 'NR > 1 { if (\$6 == "$cell_type") { print \$2; } }' $params.metadata | sort | uniq > samples.txt

        grep -E "^[1-9]|^X" $peaks > cell_type_peaks.bed
        python3 $workflow.projectDir/peak_shape_features.py cell_type_peaks.bed samples.txt ${cell_type}
        """
}

//-----------------------------------------------------
// Workflow
//-----------------------------------------------------

workflow {

    CONSENSUS_SHAPE_FEATURES()

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

    CELL_TYPE_SHAPE_FEATURES(samples_list.peaks, samples_list.cell_type)
}
