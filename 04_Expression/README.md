# Expression

This folder contains the scripts for the co-expression analysis. The order of running the scripts is described below.

1. [WGCNA](#wgcna)
    1. [WGCNA on Gene Expression](#wgcna-on-gene-expression)
    2. [Module Annotation](#module-annotation)
    3. [Module Associations with SRS and Cell Frequencies](#module-associations-with-srs-and-cell-frequencies)
    4. [Module Associations with Time Point and Diagnosis](#module-associations-with-time-point-and-diagnosis)
    5. [Module Associations with Outcome](#module-associations-with-outcome)
    6. [Module Associations](#module-associations)
    7. [Network Analysis](#network-analysis)
    8. [Preparation for Module QTL Mapping](#preparation-for-module-qtl-mapping)
    9. [Prepare Genotypes for Module QTL](#prepare-genotypes-for-module-qtl)
    10. [Perform Module QTL Mapping](#perform-module-qtl-mapping)
    11. [Module QTL Analysis](#module-qtl-analysis)
    12. [Module QTL Information](#module-qtl-information)
    13. [Extract Genotypes for Module QTL Loci Summary Statistics](#extract-genotypes-for-module-qtl-loci-summary-statistics)
    14. [Generate Module QTL Loci Summary Statistics](#generate-module-qtl-loci-summary-statistics)
    15. [Module QTL Replication](#module-qtl-replication)
    16. [Prediction from Module Eigengenes](#prediction-from-module-eigengenes)
2. [Trajectories](#trajectories)
    1. [DDRTree](#ddrtree)
    2. [Pseudotime Analysis](#pseudotime-analysis)

---

## WGCNA

For the WGCNA analysis, run the following scripts in order:

---

### WGCNA on Gene Expression

**File**: `wgcna_01_gene_expression.ipynb`

**Description**: Use the WGCNA R package on the gene expression matrix to identify modules.

**Expected Inputs**:
1. logCPM Gene Expression Matrix
2. RNA-seq Sample Information
3. Gene Information
4. Sample Key

**Dependencies**: None

**Outputs**:
1. Correlation Matrix
2. Adjacency Matrix
3. TOM Matrix
4. Gene Assignments to Modules
5. Eigengenes
6. Network Connectivity Metrics

---

### Module Annotation

**File**: `wgcna_02_gene_coexpression_module_annotation.ipynb`

**Description**: Perform module annotation using pathway information from GO, KEGG, and Reactome.

**Expected Inputs**:
1. Outputs from `wgcna_01_gene_expression.ipynb`
2. Gene Information

**Dependencies**: None

**Outputs**:
1. GO Cellular Components, Molecular Functions, and Biological Processes
2. KEGG Human Pathways
3. Reactome Pathways

---

### Module Associations with SRS and Cell Frequencies

**File**: `wgcna_03_gene_coexpression_module_srs_and_xcell.ipynb`

**Description**: Tests for association of module eigengenes with SRS assignment, cell proportions measured in the cohort, and cell frequencies derived using CIBERSORTx with the Sepsis Immunomics reference scRNA-seq data.

**Expected Inputs**:
1. Outputs from `wgcna_01_gene_expression.ipynb`
2. logCPM Gene Expression Matrix
3. Gene Information
4. Sample Information
5. SRS Assignments for Samples
6. Cell Proportions for Samples
7. xCell Information
8. CIBERSORTx Inferred Frequencies
9. scRNA-seq Cell Types
10. scRNA-seq Gene Markers

**Dependencies**: None

**Outputs**:
1. SRSq Associations with Module Eigengenes (Spearman's Rho)
2. Cell Proportion Associations with Module Eigengenes (Spearman's Rho)
3. Cell Frequency Associations with Module Eigengenes (Spearman's Rho)
4. xCell Enrichment Score Associations with Module Eigengenes (Spearman's Rho)
5. xCell Signature Enrichment in Module Eigengenes (Hypergeometric Test)
6. scRNA-seq Marker Gene Set Enrichment in Module Eigengenes (Hypergeometric Test)

---

### Module Associations with Time Point and Diagnosis

**File**: `wgcna_04_gene_coexpression_module_time_point_and_diagnosis.ipynb`

**Description**: Tests for association of module eigengenes with time point (D1, D3, D5) and diagnosis (CAP vs. FP).

**Expected Inputs**:
1. Outputs from `wgcna_01_gene_expression.ipynb`
2. Sample Information
3. Clinical Data for Patients

**Dependencies**: None

**Outputs**:
1. Time Point Association with Module Eigengenes (Repeat Measures ANOVA - GES)
2. Diagnosis Association with Module Eigengenes (Repeat Measures ANOVA - GES)

---

### Module Associations with Outcome

**File**: `wgcna_05_gene_coexpression_module_proportional_hazards_model.ipynb`

**Description**: Tests for association of module eigengenes with 28-day outcome.

**Expected Inputs**:
1. Outputs from `wgcna_01_gene_expression.ipynb`
2. Outcome Data for Patients
3. Sample Information

**Dependencies**: None

**Outputs**:
1. Outcome Association with Module Eigengenes (Cox Proportional Hazards Model - Beta)

---

### Module Associations

**File**: `wgcna_06_gene_coexpression_module_associations.ipynb`

**Description**: Compiles all association results for the module eigengenes.

**Expected Inputs**:
1. Outputs from `wgcna_01_gene_expression.ipynb`
2. Outputs from `wgcna_03_gene_coexpression_module_srs_and_xcell.ipynb`
3. Outputs from `wgcna_04_gene_coexpression_module_time_point_and_diagnosis.ipynb`
4. Outputs from `wgcna_05_gene_coexpression_module_proportional_hazards_model.ipynb`
5. xCell Information
6. scRNA-seq Cell Types

**Dependencies**: None

**Outputs**:
1. Combined Association Table

---

### Network Analysis

**File**: `wgcna_07_network_analysis.ipynb`

**Description**: Identify hub genes and export module networks for analysis with Cytoscape.

**Expected Inputs**:
1. Outputs from `wgcna_01_gene_expression.ipynb`
2. Gene Information

**Dependencies**: None

**Outputs**:
1. Module Adjacencies
2. Module Node Lists

---

### Preparation for Module QTL Mapping

**File**: `wgcna_08_eigengene_sva_preparation.ipynb`

**Description**: Prepare files for the module QTL mapping. Required before running the Nextflow pipelines.

**Expected Inputs**:
1. Outputs from `wgcna_01_gene_expression.ipynb`
2. Clinical Data for Patients
3. SRS Assignment for Samples
4. eQTL Mapping Covariates (Genotyping PCs, PEER Factors, Cell Proportions, Diagnosis, SRS Status)
5. logCPM Gene Expression Matrix
6. cis-eQTL
7. Conditional cis-eQTL
8. EBI GWAS Catalog Studies
9. EBI GWAS Catalog Associations
10. Genotyped Variants

**Dependencies**: None

**Outputs**:
1. Patient IDs for Mapping
2. List of Eigengenes
3. Mapping Data (Eigengenes, Covariates)
4. List of SNPs to Test

---

### Prepare Genotypes for Module QTL

**File**: `extract_genotypes.sh`

**Description**: Extract the genotypes for the ~70,000 SNPs that are tested for module QTL.

**Expected Inputs**:
1. Outputs from `wgcna_08_eigengene_sva_preparation.ipynb`
2. Patient Genotypes

**Dependencies**:
1. `extract_genotypes.nf`

**Outputs**:
1. Extracted, filtered, recoded genotypes for mapping

---

### Perform Module QTL Mapping

**File**: `eigengene_single_variant_association.sh`

**Description**: Perform module QTL mapping by running the linear mixed model for all the tested SNPs.

**Expected Inputs**:
1. Outputs from `wgcna_08_eigengene_sva_preparation.ipynb`
2. Outputs from `extract_genotypes.sh`

**Dependencies**:
1. `eigengene_single_variant_association.nf`
2. `eigengene_single_variant_association.R`

**Outputs**:
1. Summary Statistics for Tested SNPs for each Module Eigengene

---

### Module QTL Analysis

**File**: `wgcna_09_eigengene_sva_analysis.ipynb`

**Description**: Analyse the module QTL that are detected in the tested SNPs.

**Expected Inputs**:
1. Outputs from `wgcna_01_gene_expression.ipynb`
2. Outputs from `wgcna_08_eigengene_sva_preparation.ipynb`
3. Outputs from `eigengene_single_variant_association.sh`
4. Patient Genotyped Variants
4. Gene Information

**Dependencies**: None

**Outputs**:
1. Lead Module QTL
2. All Module QTL
3. Module QTL Full Summary Statistics
4. mQTL RDS Object

---

### Module QTL Information

**File** `wgcna_10_eigengene_sva_qtl_info.ipynb`

**Description**: Describe the proportion of input SNPs that were module QTL and identify traits associated with module QTL.

**Expected Inputs**:
1. Outputs from `wgcna_01_gene_expression.ipynb`
2. Outputs from `wgcna_08_eigengene_sva_preparation.ipynb`
3. Outputs from `wgcna_09_eigengene_sva_analysis.ipynb`
4. Gene Information
5. EBI GWAS Studies

**Dependencies**: None

**Outputs**:
1. EBI GWAS Variants that are Module QTL

---

### Extract Genotypes for Module QTL Loci Summary Statistics

**File**: `extract_ss_genotypes.sh`

**Description**: Extracts genotypes around module QTL to generate full summary statistics in the region. Also extracts genotypes for the module QTL replication.

**Expected Inputs**:
1. Outputs from `wgcna_08_eigengene_sva_preparation.ipynb`
2. Outputs from `wgcna_09_eigengene_sva_analysis.ipynb`
3. Patient Genotypes

**Dependencies**:
1. `extract_ss_genotypes.nf`

**Outputs**:
1. Extracted, filtered, recoded genotypes for module QTL loci
2. Extracted, filtered, recoded genotypes for module QTL replication

---

### Generate Module QTL Loci Summary Statistics

**File**: `module_qtl_summary_statistics.sh`

**Description**: Generates summary statistics around the detected module QTL locus.

**Expected Inputs**:
1. Outputs from `wgcna_08_eigengene_sva_preparation.ipynb`
2. Outputs from `wgcna_09_eigengene_sva_analysis.ipynb`
3. Outputs from `extract_ss_genotypes.sh`

**Dependencies**:
1. `module_qtl_summary_statistics.nf`
2. `module_qtl_summary_statistics.R`
3. `module_qtl_summary_statistics_all_pcs.R`

**Outputs**:
1. Summary statistics for module eigengenes from module QTL loci

---

### Module QTL Replication

**File**: `wgcna_11_module_qtl_replication.ipynb`

**Description**: Test the module QTL for replication in the GAinS microarray data.

**Expected Inputs**:
1. Outputs from `wgcna_01_gene_expression.ipynb`
2. Outputs from `wgcna_09_eigengene_sva_analysis.ipynb`
3. Outputs from `extract_ss_genotypes.sh`
3. Microarray Expression Data

**Dependencies**: None

**Outputs**: None

---

### Prediction from Module Eigengenes

> **Note**: Experimental

**File**: `wgcna_12_predictors.ipynb`

**Description**: Use an SVM with the module eigengenes to try and predict 28-day outcome.

**Expected Inputs**:
1. Outputs from `wgcna_01_gene_expression.ipynb`
2. Clinical Data for Patients
3. Outcome Data for Patients

**Dependencies**: None

**Outputs**: None

---

## Trajectories

For the trajectories analysis, run the following scripts in order:

---

### DDRTree

**File**: `trajectories_01_ddrtree.ipynb`

**Description**: Runs the DDRTree algorithm on the bulk RNA-seq expression data to embed patients along a trajectory and calculate pseudotime. This method was originally developed for scRNA-seq data.

**Expected Inputs**:
1. logCPM Gene Expression Matrix
2. Module Eigengenes
3. RNA-seq Sample Information
4. Module Eigengene Associations

**Dependencies**:
1. `04_Expression/scripts/utils/tree_dfs.cpp`
2. `04_Expression/scripts/utils/tree_projection.cpp`
3. `04_Expression/scripts/utils/tree_pseudotime.cpp`

**Outputs**:
1. DDRTree outputs for: (1) full expression matrix, (2) eigengene matrix, (3) eigengenes associated with time point only
2. DDRTree branches
3. DDRTree projected points
4. Pseudotime

---

### Pseudotime Analysis

**File**: `trajectories_02_pseudotime.ipynb`

**Description**: Analyses to check if pseudotime aligns patients along a clinical trajectory.

**Expected Inputs**:
1. Outputs from `trajectories_01_ddrtree.ipynb`
2. logCPM Gene Expression Matrix
3. Module Eigengenes
4. RNA-seq Sample Information
5. Clinical Data for Patients
6. SRS Assignments for Samples
7. Cell Proportions for Samples

**Dependencies**:
1. `04_Expression/scripts/utils/stree_bfs.cpp`

**Outputs**: None
