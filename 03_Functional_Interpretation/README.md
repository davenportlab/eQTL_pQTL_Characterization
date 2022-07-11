# Functional Interpretation

This folder contains the scripts for functional interpretation of the molecular QTL. The order of running the scripts is described below:

## Retrieve Public Data

Run the scripts in the following order:

### Generate Metadata

**File**: `setup_01_process_metadata.ipynb`

**Description**: Takes metadata files from the SRA and generates tables that will be used for downstream analysis.

**Expected Inputs**:
1. SRA Metadata from Studies

**Dependencies**: None

**Outputs**:
1. RNA-seq Reads Table
2. ATAC-seq Reads Table

### RNA-seq Data Overview

**File**: `rna_01_rna_seq_raw_data.ipynb`

**Description**: Overview of metadata for RNA-seq.

**Expected Inputs**:
1. Outputs from `setup_01_process_metadata.ipynb`

**Dependencies**: None

**Outputs**: None

### ATAC-seq Data Overview

**File**: `atac_01_atac_seq_raw_data.ipynb`

**Description**: Overview of metadata for ATAC-seq.

**Expected Inputs**:
1. Outputs from `setup_01_process_metadata.ipynb`

**Dependencies**: None

**Outputs**: None

### Download from SRA

**File**: `download_from_sra.sh`

**Description**: Download FASTQ files from SRA for relevant studies.

**Expected Inputs**:
1. Outputs from `setup_01_process_metadata.ipynb`

**Dependencies**:
1. `download_from_sra.nf`

**Outputs**:
1. FASTQ files from RNA-seq studies
2. FASTQ files from ATAC-seq studies

## Alignment

Run the scripts in the following order:

### Index Genomes

**File**: `index_genomes.sh`

**Description**: Index genomes for use with STAR and Bowtie 2.

**Expexted Inputs**:
1. Ensembl Genome Assembly
2. Ensembl Genome Annotation

**Dependencies**:
1. `index_genome.nf`

**Outputs**:
1. STAR Index
2. Bowtie 2 Index

## RNA Sequencing

Run the scripts in the following order:

### Trim Adapter Contamination

**File**: `01_pre_process_paired_rna_seq.sh`

**Description**: Run TrimGalore to remove adapter sequences.

**Expexcted Inputs**
1. Outputs from `download_from_sra.sh`

**Dependencies**:
1. `01_pre_process_paired_rna_seq.nf`

**Outputs**:
1. Trimmed FASTQ Files
2. MultiQC Report After Trimming

### Alignment

**File**: `02_process_paired_rna_seq.sh`

**Description**: Align reads to the genome using STAR and perform post-alignment quality control.

**Expected Inputs**:
1. Outputs from `index_genomes.sh`
2. Outputs from `01_pre_process_paired_rna_seq.sh`

**Dependencies**:
1. `02_process_paired_rna_seq.nf`

**Outputs**:
1. Alignment BAM Files
2. Read Counts in Genes from STAR
3. MultiQC Report After Post-Alignment Quality Control

### Merge Replicates and Count Reads

**File**: `03_merge_replicates.sh`

**Description**: Merge any technical replicates and then count the number of reads that overlap gene features.

**Expected Inputs**:
1. Outputs from `setup_01_process_metadata.ipynb`
2. Outputs from `02_process_paired_rna_seq.sh`
3. Ensembl Genome Annotation

**Dependencies**:
1. `03_merge_replicates.nf`

**Outputs**:
1. Merged BAM Files
2. Reads in Exons from featureCounts

### RNA-seq Analysis

**File**: `04_rna_seq_analysis.sh`

**Description**: Create an aggregate count matrix across samples and run RNASeQC.

**Expected Inputs**:
1. Outputs from `setup_01_process_metadata.ipynb`
2. Outputs from `03_merge_replicates.sh`
3. Collapsed Ensembl Genome Annotation

**Dependencies**:
1. `04_rna_seq_analysis.nf`

**Outputs**:
1. Read Count Matrix
2. Aggregate RNASeQC Output

### Basic Alignment QC

**File**: `rna_02_rna_seq_read_counts.ipynb`

**Description**: Analyse basic RNA-seq quality control metrics.

**Expected Inputs**:
1. Outputs from `setup_01_process_metadata.ipynb`
2. Outputs from `04_rna_seq_analysis.sh`

**Dependencies**: None

**Outputs**: None

### Differentially Expressed Genes

**File**: `rna_03_rna_seq_count_normalization.ipynb`

**Description**: Identify differentially expressed genes from the count matrix.

**Expected Inputs**:
1. Outputs from `setup_01_process_metadata.ipynb`
2. Outputs from `04_rna_seq_analysis.sh`
3. Ensembl Genome Annotation

**Dependencies**: None

**Outputs**:
1. Differentially Expressed Genes

## ATAC Sequencing

Run the scripts in the following order:

### Trim Adapter Contamination

**File**: `01_pre_process_paired_atac_seq.sh`

**Description**: Run TrimGalore to remove adapter sequences.

**Expexcted Inputs**
1. Outputs from `download_from_sra.sh`

**Dependencies**:
1. `01_pre_process_paired_atac_seq.nf`

**Outputs**:
1. Trimmed FASTQ Files
2. MultiQC Report After Trimming

### Alignment

**File**: `02_process_paired_atac_seq.sh`

**Description**: Align reads to the genome using Bowtie 2 and perform post-alignment quality control.

**Expected Inputs**:
1. Outputs from `index_genomes.sh`
2. Outputs from `01_pre_process_paired_atac_seq.sh`
3. Ensembl Genome Sequence
4. ENCODE Black Listed Regions

**Dependencies**:
1. `02_process_paired_atac_seq.nf`

**Outputs**:
1. Alignment BAM Files
2. MultiQC Report After Post-Alignment Quality Control

### Merge Replicates and Call Peaks

**File**: `03_call_atac_seq_peaks.sh`

**Description**: Merge any technical replicates and then call peaks using MACS2.

**Expected Inputs**:
1. Outputs from `setup_01_process_metadata.ipynb`
2. Outputs from `index_genomes.sh`
3. Outputs from `02_process_paired_atac_seq.sh`

**Dependencies**:
1. `03_call_atac_seq_peaks.nf`

**Outputs**:
1. Merged BAM Files
2. Peaks from MACS2

### Quality Control

**File**: `04_sample_quality.sh`

**Description**: Performs peak quality control by measuring peak widths and measuring Jaccard similarity between accessibility profiles. Generates TSS enrichment scores for the samples.

**Expected Inputs**:
1. Outputs from `setup_01_process_metadata.ipynb`
2. Outputs from `index_genomes.sh`
3. Outputs from `03_call_atac_seq_peaks.sh`
4. Ensembl Genome Annotation

**Dependencies**:
1. `04_sample_quality.nf`

**Outputs**:
1. Peak Widths
2. Jaccard Values between Samples
3. Sample TSS Enrichment Scores

### Peak Sets

**File**: `05_consensus_peak_set.sh`

**Description**: Generate group, cell type, and consensus peak sets from the samples.

**Expected Inputs**:
1. Outputs from `setup_01_process_metadata.ipynb`
2. Outputs from `03_call_atac_seq_peaks.sh`

**Dependencies**:
1. `05_consensus_peak_set.nf`

**Outputs**:
1. Group Peak Sets
2. Cell Type Peak Sets
3. Consensus Peak Sets
4. Group Jaccard Values
5. Cell Type Jaccard Values

### Count Fragments in Peaks

**File**: `06_count_fragments.sh`

**Description**: Count fragments in all the peak sets.

**Expected Inputs**:
1. Outputs from `setup_01_process_metadata.ipynb`
2. Outputs from `05_consensus_peak_set.sh`

**Dependencies**:
1. `06_count_fragments.nf`

**Outputs**:
1. Group Peak Sets Fragment Counts
2. Cell Type Peak Sets Fragment Counts
3. Consensus Peak Sets Fragment Counts
4. Fraction of Reads in Peaks (FRiP)

### Sample Jaccard Similarity

**File**: `atac_02_peaks_jaccard_similarity.ipynb`

**Description**: Compare peaks from individual samples using the Jaccard metric and performs MDS.

**Expected Inputs**:
1. Outputs from `setup_01_process_metadata.ipynb`
2. Outputs from `04_sample_quality.sh`

**Dependencies**: None

**Outputs**: None

### Explore Consensus Peak Set

**File**: `atac_03_consensus_peak_set.ipynb`

**Description**: Explore the peaks in the consensus peak set.

**Expected Inputs**:
1. Outputs from `setup_01_process_metadata.ipynb`
2. Outputs from `index_genomes.sh`
3. Outputs from `05_consensus_peak_set.sh`

**Dependencies**: None

**Outputs**: None

### Explore All Peak Set

**File**: `atac_04_all_peak_sets.ipynb`

**Description**: Explore the peaks in the group and cell type peak sets and compare the sets using the Jaccard metric.

**Expected Inputs**:
1. Outputs from `05_consensus_peak_set.sh`
2. Calderon et al. Lineages

**Dependencies**: None

**Outputs**: None

### Sample Alignment QC

**File**: `atac_05_frip_and_tss_enrichment.ipynb`

**Description**: Explore the fraction of reads in peaks (FRiP) and TSS enrichment scores for samples based on division by various attributes.

**Expected Inputs**:
1. Outputs from `setup_01_process_metadata.ipynb`
2. Outputs from `04_sample_quality.sh`
3. Outputs from `06_count_fragments.sh`

**Dependencies**: None

**Outputs**: None

### Differentially Accessible Peaks in Immune Atlas

**File**: `atac_06_immune_peak_count_normalization.ipynb`

**Description**: Identify differentially accessible (DA) peaks in the immune atlas.

**Expected Inputs**:
1. Outputs from `setup_01_process_metadata.ipynb`
2. Outputs from `04_sample_quality.sh`
3. Outputs from `06_count_fragments.sh`

**Dependencies**: None

**Outputs**:
1. Differentially Accessible Peaks in Immune Atlas

### Differentially Accessible Peaks in Neutrophil Atlas

**File**: `atac_07_neutrophil_peak_count_normalization.ipynb`

**Description**: Identify differentially accessible (DA) peaks in the neutrophil atlas.

**Expected Inputs**:
1. Outputs from `setup_01_process_metadata.ipynb`
2. Outputs from `04_sample_quality.sh`
3. Outputs from `06_count_fragments.sh`

**Dependencies**: None

**Outputs**:
1. Differentially Accessible Peaks in Neutrophil Atlas

### Compare with Original Calderon et al. Analysis

**File**: `atac_08_calderon_et_al_comparison.ipynb`

**Description**: Compare the analysis conducted with the original Calderon et al. fragment counts and DA peaks.

**Expected Inputs**:
1. Outputs from `setup_01_process_metadata.ipynb`
2. Outputs from `06_count_fragments.sh`
3. Outputs from `atac_06_immune_peak_count_normalization.ipynb`
4. Original Calderon et al. Peak Count Matrix
5. Original Calderon et al. Differentially Accessible Peaks
6. Chain from hg19 to GRCh38

**Dependencies**: None

**Outputs**: None

### Motif Enrichment

**File**: `07_motifs.sh`

**Description**: Run motif enrichment and peak annotation analysis on peak sets and DA peaks.

**Expected Inputs**:
1. Outputs from `index_genomes.sh`
2. Outputs from `05_consensus_peak_set.sh`
3. Outputs from `atac_06_immune_peak_count_normalization.ipynb`
4. Outputs from `atac_07_neutrophil_peak_count_normalization.ipynb`
5. Ensembl Genome Annotation
6. Ensembl Genome Sequence
7. JASPAR Motifs

**Dependencies**:
1. `07_motifs.nf`

**Outputs**:
1. HOMER Peak Annotations for Consensus Peak Sets
2. HOMER Peak Annotations for Cell Type Peak Sets
3. HOMER Peak Annotations for Group Peak Sets
4. Motif Enrichment for Consensus Peak Sets
5. Motif Enrichment for Cell Type Peak Sets
6. Motif Enrichment for Group Peak Sets
7. Motif Enrichment for DA Peak Sets (More and Less Accessible)

### Analysis of Motif Enrichment

**File**: `atac_09_homer.ipynb`

**Description**: Analyse results from the motif enrichment and peak annotation.

**Expected Inputs**:
1. Outputs from `07_motifs.sh`

**Dependencies**: None

**Outputs**: None

## Functional Enrichment

Run the scripts in the following order:

### Gene-Peak Clusters

> **Note**: Experimental

**File**: `enrichment_08_gene_peak_clusters.ipynb`

**Description**: Identify modules of accessibility by clustering peaks near genes based on the count data.

**Expected Inputs**:
1. Outputs from `setup_01_process_metadata.ipynb`
2. Outputs from `06_count_fragments.sh`
3. Ensembl Genome Annotation

**Dependencies**: None

**Outputs**: None

### Overlap with Prior Analyses

> **Note**: Experimental

**File**: `enrichment_07_overlap_with_prior_analyses.ipynb`

**Description**: Overlap of cis-eQTL with results from Calderon et al. and Ram-Mohan et al.

**Expected Inputs**:
1. Calderon et al. Differentially Accessible Peaks
2. Calderon et al. Chromatin Accessibility QTL (caQTL)
3. Ram-Mohan et al. Differentially Accessible Peaks
4. Cis-eQTL
5. Conditional cis-eQTL
6. Patients Genotyped Variants
7. Gene Information
8. Chain from GRCh38 to hg19

**Dependencies**: None

**Outputs**:
1. Cis-eQTL and Tagging SNPs on hg19
2. Conditional cis-eQTL and Tagging SNPs on hg19
3. TSS Regions on hg19
4. SNPs Overlapping Differentially Accessible Peaks
5. SNPs that are caQTL

### SNP Matching Within Cohort

> **Note**: Experimental

**File**: `snp_matching.sh`

**Description**: This is an attempt to perform SNP matching using LD information from within the cohort.

**Expected Inputs**:
1. Outputs from `ld_fine_mapping.sh`
2. Outputs from `05_consensus_peak_set.sh`
3. SNPsnap Results
4. Gene Information
5. Patient Genotyped Variants
6. ENCODE cCREs
7. ChromHMM States
8. Sepsis-Enhanced eSNPs

**Dependencies**:
1. `snp_matching.nf`
2. `snp_matching_generate_intervals.R`
3. `snp_matching_tests.R`

**Outputs**:
1. Fisher Test Results from Enrichment
2. TSS-Corrected Enrichment Test Results

### Analysis of SNP Matching Within Cohort

> **Note**: Experimental

**File**: `enrichment_10_snp_matching.ipynb`

**Description**: Analyse results of running SNP matching within the cohort.

**Expected Inputs**:
1. Outputs from `snp_matching.sh`

**Dependencies**: None

**Outputs**: None

### SNP Matching Setup

**File**: `enrichment_01_snp_snap_setup.ipynb`

**Description**: Export SNPs to be used as input for the SNPsnap server.

**Expected Inputs**:
1. Outputs from `ld_fine_mapping.sh`
1. Cis-eQTL
2. Conditional cis-eQTL

**Dependencies**: None

**Outputs**
1. Table of All Cis-eQTL eSNPs

### Matched SNP Enrichment

**File**: `snp_snap_matching.sh`

**Description**: Measures the number of peaks that overlap SNPs from various SNP lists.

**Expected Inputs**:
1. Outputs from `ld_fine_mapping.sh`
2. Outputs from `05_consensus_peak_set.sh`
3. SNPsnap Results
4. Patient Genotyped Variants
5. ENCODE cCREs
6. ChromHMM States
7. Sepsis-Enhanced eSNPs
8. Chain from hg19 to GRCh38

**Dependencies**:
1. `snp_snap_matching.nf`

**Outputs**:
1. Cis-eQTL Overlaps
2. Conditional cis-eQTL Overlaps
3. Sepsis-Enhanced eQTL Overlaps
4. Background Overlaps

### Analysis of Matched SNP Enrichment

**File**: `enrichment_02_snp_snap.ipynb`

**Description**: Analysis of enrichment results from the matched SNP approach.

**Expected Inputs**:
1. Outputs from `snp_snap_matching.sh`
2. Calderon et al. Lineages
3. Roadmap Epigenomics ChromHMM States and Epigenomes

**Dependencies**: None

**Outputs**: None

### CHEERS on Calderon et al. Data

> **Note**: Experimental

**File**: `cheers_calderon_et_al_hg19.sh`

**Description**: Perform CHEERS on peaks generated by the Calderon et al. analysis.

**Expected Inputs**:
1. Outputs from `ld_fine_mapping.sh`
2. Outputs from `enrichment_07_overlap_with_prior_analyses.ipynb`
3. Calderon et al. Peak Count Matrix
4. eQTL Data from Fairfax et al. 
5. Chain from hg18 to hg19
6. Conditional cis-eQTL
7. Sepsis-Enhanced eQTL

**Dependencies**:
1. `cheers_calderon_et_al_hg19.nf`
2. `calderon_et_al_hg19_positive_controls.R`

**Outputs**:
1. Normalised Peak Specificity Scores
2. CHEERS Results

### Analysis of CHEERS on Calderon et al. Data

> **Note**: Experimental

**File**: `enrichment_06_cheers_hg19_calderon_et_al.ipynb`

**Description**: Analyse results from CHEERS on the Calderon et al. data.

**Expected Inputs**:
1. Outputs from `cheers_calderon_et_al_hg19.sh`
2. Calderon et al. Lineages

**Dependencies**: None

**Outputs**: None

### CHEERS

**File**: `cheers.sh`

**Description**: Run CHEERS on various SNP lists in the neutrophil atlas samples.

**Expected Inputs**:
1. Outputs from `ld_fine_mapping.sh`
2. Outputs from `06_count_fragments.sh`
3. Patients Genotyped Variants
4. Gene Information
5. Conditional cis-eQTL
6. Sepsis-Enhanced eQTL

**Dependencies**:
1. `cheers.nf`

**Outputs**:
1. Normalised Peak Specificity Scores
2. CHEERS Results for Ex Vivo Stimulations
3. CHEERS Results for Whole Blood Stimulations

### Analysis of CHEERS

**File**: `enrichment_03_cheers.ipynb`

**Description**: Analyse results from CHEERS.

**Expected Inputs**:
1. Outputs from `cheers.sh`

**Dependencies**: None

**Outputs**: None

### GoShifter on Calderon et al. Data

> **Note**: Experimental

**File**: `goshifter_calderon_et_al_hg19.sh`

**Description**: Perform GoShifter on peaks from original Calderon et al. analysis.

**Expected Inputs**:
1. Outputs from `ld_fine_mapping.sh`
2. Outputs from `eqtl_in_epigenetic_data.ipynb`
3. Calderon et al. Peak Count Matrix
4. Calderon et al. Differentially Accessible Peaks
5. eQTL Data from Fairfax et al. 
6. Chain from hg18 to hg19
7. Conditional cis-eQTL
8. Sepsis-Enhanced eQTL

**Dependencies**:
1. `goshifter_calderon_et_al_hg19.nf`

**Outputs**:
1. Overlap Scores
2. Empirical P-Values

### Analysis of GoShifter on Calderon et al. Data

> **Note**: Experimental

**File**: `enrichment_09_goshifter_hg19_calderon_et_al.ipynb`

**Description**: Analyse results from GoShifter on the Calderon et al. data.

**Expected Inputs**:
1. Outputs from `goshifter_calderon_et_al_hg19.sh`
2. Calderon et al. Lineages

**Dependencies**: None

**Outputs**: None

### GoShifter

**File**: `go_shifter.sh`

**Description**: Perform GoShifter analysis for various SNP lists in the accessibility peak sets.

**Expected Inputs**:
1. Outputs from `ld_fine_mapping.sh`
2. Outputs from `05_consensus_peak_set.sh`
3. Gene Information
4. Patients Genotyped Variants
5. Conditional cis-eQTL
6. Sepsis-Enhanced eQTL

**Dependencies**:
1. `go_shifter.nf`

**Outputs**:
1. SNP Lists
2. Peak Lists
3. Overlap Scores
4. Empirical P-Values

### Analysis of GoShifter

**File**: `enrichment_04_goshifter.ipynb`

**Description**: Analysis of GoShifter results.

**Expected Inputs**:
1. Outputs from `setup_01_process_metadata.ipynb`
2. Outputs from `go_shifter.sh`
3. Conditional cis-eQTL
4. Gene Information
5. Calderon et al. Lineages

**Dependencies**: None

**Outputs**: None

### Create GRMs

**File**: `heritability.sh`

**Description**: Generate genetic relationship matrices (GRMs) using GCTA based on all the peak sets.

**Expected Inputs**:
1. Outputs from `wgcna_08_eigengene_sva_preparation.ipynb`
2. Outputs from `05_consensus_peak_set.sh`
3. ENCODE cCREs
4. ChromHMM States
5. Patient Genotypes

**Dependencies**:
1. `heritability.nf`

**Outputs**:
1. Genotypes for SNPs on Autosomes
2. Full Genetic Relationship Matrix
3. Genetic Relationship Matrices for Genome Annotations

### Partitioned Heritability

**File**: `paritioned_heritability.sh`

**Description**: Run variance components model to partition heritability based on various annotations of the genome.

**Expected Inputs**:
1. Outputs from `wgcna_01_gene_expression.ipynb`
2. Outputs from `wgcna_08_eigengene_sva_preparation.ipynb`
3. Outputs from `heritability.sh`

**Dependencies**:
1. `partitioned_heritability.nf`
2. `partitioned_heritability.R`

**Outputs**:
1. Proportion of Variance Assigned to Each Component

### Analysis of Partitioned Heritability

**File**: `enrichment_05_partitioned_heritability.ipynb`

**Description**: Analyse results from the variance components model.

**Expected Inputs**:
1. Outputs from `heritability.sh`
2. Outputs from `paritioned_heritability.sh`
3. Calderon et al. Lineages

**Dependencies**: None

**Outputs**: None

### Scissor scRNA-seq Analysis

> **Note**: Experimental

**File**: `enrichment_11_scissors_scrna_seq.ipynb`

**Description**: Use Scissor to identify cell types that are important to outcome directly from scRNA-seq data in whole blood.

**Expected Inputs**:
1. COMBAT scRNA-seq Data Set
2. logCPM Gene Expression Matrix
3. RNA-seq Sample Information
4. RNA-seq Sample Key
5. RNA-seq Outcome Information

**Dependencies**: None

**Output**: None

## Variant Effect Prediction

Run the scripts in the following order:

### Preparation

**File**: `vep_01_qtl_vep_analyses.ipynb`

**Description**: Export SNP lists that will be used for variant effect prediction.

**Expected Inputs**:
1. Outputs from `wgcna_01_gene_expression.ipynb`
2. Outputs from `extract_genotypes.sh`
3. Outputs from `wgcna_09_eigengene_sva_analysis.ipynb`
4. Outputs from `extract_pqtl.sh`
5. Patients Genotyped Variants
6. Cis-eQTL
7. Conditional cis-eQTl
8. SuSiE Credible Sets
9. Chain from GRCh38 to hg19

**Dependencies**: None

**Outputs**:
1. VCF with Variants with Both Alleles as the Reference
2. VCF with QTL Variants on hg19
3. VCF with QTL Variants on hg19 with CS Variants
4. VCF with QTL Variants on hg19 with CS Variants with Both Alleles as the Reference

### Ensembl's VEP

**File**: `ensembl_vep.sh`

**Description**: Run Ensembl's VEP that is maintained by Sanger's HGI.

**Expected Inputs**:
1. Outputs from `vep_01_qtl_vep_analyses.ipynb`

**Dependencies**:
1. `ensembl_vep_post_processing.py`

**Outputs**:
1. Table of VEP Consequences

### VEP Analysis

**File**: `vep_03_ensembl_vep.ipynb`

**Description**: Analysis of results from variant effect prediction.

**Expected Inputs**:
1. Outputs from `wgcna_09_eigengene_sva_analysis.ipynb`
2. Outputs from `extract_pqtl.sh`
3. Outputs from `ensembl_vep.sh`
4. Gene Information
5. Conditional cis-eQTL

**Dependencies**: None

**Outputs**: None

## Regulation

Run the scripts in the following order:

### Peak Shapes

> **Note**: Experimental

**File**: `peak_shape_features.sh`

**Description**: Generate features that describe the shape of the feature.

**Expected Inputs**:
1. Outputs from `setup_01_process_metadata.ipynb`
2. Outputs from `05_consensus_peak_set.sh`

**Dependencies**:
1. `peak_shape_features.nf`
2. `peak_shape_features.py`

**Outputs**:
1. Fiedler Vectors for Consensus and Cell Type Peaks
2. Coverage PCs for Consensus and Cell Type Peaks

## Integration

Run the scripts in the following order:

### Module Integration

**File**: `integration_01_modules.ipynb`

**Description**: Integrate information from all the analyses for some interesting modules.

**Expected Inputs**:
1. Outputs from `wgcna_01_gene_expression.ipynb`
2. Outputs from `04_cis_eqtl_module_qtl_all_pcs_colocalization.ipynb`
3. Outputs from `ensembl_vep.sh`
4. Outputs from `go_shifter.sh`
5. Outputs from `paritioned_heritability.sh`
6. Outputs from `mqtl_susie_fine_mapping.sh`
7. Gene Information
8. Patients Genotyped Variants
9. Conditional cis-eQTL
10. Calderon et al. Lineages

**Dependencies**: None

**Outputs**: None
