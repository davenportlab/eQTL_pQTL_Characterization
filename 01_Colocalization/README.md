# Colocalisation

This folder contains the scripts for colocalisation and fine mapping. The order of running the scripts is described below:

## Conditional Summary Statistics

Run the scripts in the following order:

### Cis-eQTL Conditional Summary Statistics

**File**: `eqtl_conditional_effects.sh`

**Description**: Run the eQTL mapping LMM, conditioning on signals to generate summary statistics for each independent signal.

**Expected Inputs**:
1. Patient Genotypes
2. Gene Information
3. Conditional cis-eQTL
4. cis-eQTL Covariates (Genotyping PCs, PEER Factors, etc.)

**Dependencies**:
1. `eqtl_conditional_effects.nf`
2. `eqtl_conditional_effects_prepare_loci.R`
3. `eqtl_conditional_effects_association.R`

**Outputs**:
1. Conditional summary statistics for individual conditional cis-eQTL signals

## LD Fine Mapping

This analysis was concerned with identifying variants in LD with molecular QTL. Run the scripts in the following order:

### Retrieve Tagging SNPs

**File**: `ld_fine_mapping.sh`

**Description**: Retrieves SNPs tagging various lists of SNPs.

**Expected Inputs**:
1. cis-eQTL
2. Conditional cis-eQTL
3. cis-pQTL
4. trans-pQTL
5. Module QTL

**Dependencies**:
1. `ld_fine_mapping.nf`
2. `ld_fine_mapping_retrieve_snps.R`

**Output**:
1. Tagging SNPs in 1 Mb for QTL with $R^2 = 1$
2. Tagging SNPs in 1 Mb for QTL with $R^2 > 0.9$
3. Tagging SNPs in 1 Mb for QTL with $R^2 > 0.8$

## Fine Mapping Setup

The purpose of this section is to retrieve the genotypes for the various QTL loci to generate LD matrices. Run the scripts in the following order:

### Extract Genotypes

**File**: `extract_genotypes.sh`

**Description**: Extracts genotypes for all the SNPs tested for cis-eQTL and pQTL. Only retrieves genotypes from the patients used.

**Expected Inputs**:
1. Patient Genotypes
1. cis-eQTL results
2. cis-pQTL and trans-pQTL results

**Dependencies**:
1. `extract_genotypes.nf`
2. `extract_genotypes.R`

**Outputs**:
1. cis-pQTL Genotypes
2. trans-pQTL Genotypes
3. eQTL Genotypes

## FINEMAP

Run the scripts in the following order:

### Cis-eQTL Fine Mapping

**File**: `cis_eqtl_finemap_fine_mapping.sh`

**Description**: Use FINEMAP for the cis-eQTL results.

**Expected Inputs**:
1. Outputs from `eqtl_conditional_effects.sh`
2. Outputs from `extract_genotypes.sh`
3. Cis-eQTL Summary Statistics
4. Patient Genotyped Variants
5. logCPM Gene Expression Matrix

**Dependencies**:
1. `cis_eqtl_finemap_fine_mapping.nf`
2. `cis_eqtl_finemap_fine_mapping_split_loci.R`
3. `cis_eqtl_finemap_fine_mapping_aggregate.py`

**Outputs**:
1. FINEMAP Credible Sets for cis-eQTL
2. FINEMAP Credible Sets for Conditional cis-eQTL

### Module QTL Fine Mapping

**File**: `mqtl_finemap_fine_mapping.sh`

**Description**: Use FINEMAP for the module QTL results.

**Expected Inputs**:
1. Outputs from `module_qtl_summary_statistics.sh`
2. Outputs from `extract_ss_genotypes.sh`
3. Patient Genotyped Variants
4. logCPM Gene Expression Matrix

**Dependencies**:
1. `mqtl_finemap_fine_mapping.nf`
2. `mqtl_finemap_fine_mapping_split_loci.R`
3. `mqtl_finemap_fine_mapping_aggregate.py`

**Outputs**:
1. FINEMAP Credible Sets for Module QTL

### pQTL Fine Mapping

**File**: `pqtl_finemap_fine_mapping.sh`

**Description**: Use FINEMAP for the pQTL results.

**Expected Inputs**:
1. Outputs from `extract_pqtl.sh`
2. Outputs from `extract_genotypes.sh`
4. Patient Genotyped Variants
5. Protein Sample Information

**Dependencies**:
1. `pqtl_finemap_fine_mapping.nf`
2. `pqtl_finemap_fine_mapping_split_loci.R`
3. `pqtl_finemap_fine_mapping_aggregate.py`

**Outputs**:
1. FINEMAP Credible Sets for cis-pQTL
2. FINEMAP Credible Sets for trans-pQTL

## SuSiE

Run the scripts in the following order:

### Extract 1000G LD

> **Note**: Experimental

**File**: `extract_mqtl_ld_1000g.sh`

**Description**: Extracts genotypes from the 1000G European superpopulation at module QTL loci for fine mapping of GWAS associations.

**Expected Inputs**:
1. Outputs from `module_qtl_summary_statistics.sh`
2. VCFs from 1000 Genomes Project
3. European Samples from 1000 Genomes Project

**Dependencies**:
1. `extract_mqtl_ld_1000g.nf`

**Outputs**:
1. Extracted, filtered, recoded genotypes from 1000G

### Cis-eQTL Fine Mapping

**File**: `cis_eqtl_susie_fine_mapping.sh`

**Description**: Use SuSiE for the cis-eQTL results.

**Expected Inputs**:
1. Outputs from `eqtl_conditional_effects.sh`
2. Outputs from `extract_genotypes.sh`
3. Cis-eQTL Summary Statistics
4. Patient Genotyped Variants
5. logCPM Gene Expression Matrix

**Dependencies**:
1. `cis_eqtl_susie_fine_mapping.nf`
2. `cis_eqtl_susie_fine_mapping_split_loci.R`
3. `cis_eqtl_susie_fine_mapping.R`

**Outputs**:
1. SuSiE Credible Sets for cis-eQTL
2. SuSiE Credible Sets for Conditional cis-eQTL

### Module QTL Fine Mapping

**File**: `mqtl_susie_fine_mapping.sh`

**Description**: Use SuSiE for the module QTL results.

**Expected Inputs**:
1. Outputs from `module_qtl_summary_statistics.sh`
2. Outputs from `extract_ss_genotypes.sh`
3. Outputs from `wgcna_01_gene_expression.ipynb`
4. Patient Genotyped Variants

**Dependencies**:
1. `mqtl_susie_fine_mapping.nf`
2. `mqtl_susie_fine_mapping_split_loci.R`
3. `mqtl_susie_fine_mapping.R`

**Outputs**:
1. SuSiE Credible Sets for Module QTL

### pQTL Fine Mapping

**File**: `pqtl_susie_fine_mapping.sh`

**Description**: Use SuSiE for the pQTL results.

**Expected Inputs**:
1. Outputs from `extract_pqtl.sh`
2. Outputs from `extract_genotypes.sh`
3. Patient Genotyped Variants
4. Protein Sample Information
5. Protein Expression

**Dependencies**:
1. `pqtl_susie_fine_mapping.nf`
2. `pqtl_susie_fine_mapping_split_loci.R`
3. `pqtl_susie_fine_mapping.R`

**Outputs**:
1. SuSiE Credible Sets for cis-pQTL
2. SuSiE Credible Sets for trans-pQTL

## Colocalisation

Run the scripts in the following order:

### Cis-eQTL Colocalisation

**File**: `01_cis_eqtl_colocalization.ipynb`

**Description**: Colocalise cis-eQTL with each other based on the observation that many conditional eGenes share the same lead eSNP. Generates some example figures of colocalising cis-eQTL and also checks which module QTL are present in the set of colocalising cis-eQTL.

**Expected Inputs**:
1. Outputs from `eqtl_conditional_effects.sh`
2. Outputs from `extract_genotypes.sh`
3. Outputs from `wgcna_09_eigengene_sva_analysis.ipynb`
4. Conditional cis-eQTL
5. logCPM Gene Expression Matrix
6. Gene Information

**Dependencies**: None

**Outputs**:
1. Colocalising cis-eQTL

### Cis-eQTL Colocalisation with Cis-pQTL

**File**: `02_cis_eqtl_cis_pqtl_colocalization.ipynb`

**Description**: Colocalise cis-eQTL with cognate cis-pQTL loci. Generates some example figures of colocalising loci.

**Expected Inputs**:
1. Outputs from `protein_metadata.ipynb`
2. Outputs from `extract_pqtl.sh`
3. Outputs from `eqtl_conditional_effects.sh`
4. Outputs from `extract_genotypes.sh`
5. logCPM Gene Expression Matrix
6. Protein Expression Matrix
7. Gene Information

**Dependencies**: None

**Outputs**: None

### Trans-pQTL Colocalisation

**File**: `03_trans_pqtl_colocalization.ipynb`

**Description**: Colocalise trans-pQTL with each other and with cis-eQTL and cis-pQTL. Generates figures of the colocalising loci.

**Expected Inputs**:
1. Outputs from `extract_pqtl.sh`
2. Outputs from `extract_genotypes.sh`
3. Outputs from `eqtl_conditional_effects.sh`
4. Outputs from `eqtl_conditional_effects.sh`
5. logCPM Gene Expression Matrix
6. Protein Expression Matrix

**Dependencies**: None

**Outputs**:
1. Colocalising cis-eQTL with Module QTL

**Dependencies**: None

**Outputs**: None

### Cis-eQTL Colocalising with Module QTL

**File**: `04_cis_eqtl_module_qtl_all_pcs_colocalization.ipynb`

**Description**: Colocalise cis-eQTL with module QTL.

**Expected Inputs**:
1. Outputs from `module_qtl_summary_statistics.sh`
2. Outputs from `wgcna_01_gene_expression.ipynb`
3. Outputs from `extract_genotypes.sh`
4. Patient Genotyped Variants
5. logCPM Gene Expression Matrix

**Dependencies**: None

**Outputs**: None

### Module QTL Colocalising with GWAS Associations

**File**: `05_module_qtl_gwas_colocalization.ipynb`

**Description**: Colocalise module QTL with GWAS associations from the EBI GWAS Catalog.

**Expected Inputs**:
1. Outputs from `module_qtl_summary_statistics.sh`
2. Outputs from `extract_genotypes.sh`
3. Outputs from `wgcna_01_gene_expression.ipynb`
4. Outputs from `wgcna_10_eigengene_sva_qtl_info.ipynb`
5. Patient Genotyped Variants
6. EBI GWAS Catalog Full Summary Statistics

**Dependencies**: None

**Outputs**: None

### Compare Credible Sets from Cis-eQTL

**File**: `06_compare_credible_sets.ipynb`

**Description**: Compare FINEMAP CSs, SuSiE CSs, and tagging SNP sets.

**Expected Inputs**:
1. Outputs from `ld_fine_mapping.sh`
2. Outputs from `cis_eqtl_finemap_fine_mapping.sh`
3. Outputs from `cis_eqtl_susie_fine_mapping.sh`
4. Cis-eQTL

**Dependencies**: None

**Outputs**: None

### Compare Credible Sets from Module QTL and pQTL

**File**: `07_mqtl_pqtl_credible_sets.ipynb`

**Description**: Compare FINEMAP CSs, SuSiE CSs, and tagging SNPs from module QTL and pQTL.

**Expected Inputs**:
1. Outputs from `wgcna_09_eigengene_sva_analysis.ipynb`
2. Outputs from `ld_fine_mapping.sh`
3. Outputs from `extract_pqtl.sh`
4. Outputs from `mqtl_finemap_fine_mapping.sh`
5. Outputs from `pqtl_finemap_fine_mapping.sh`
6. Outputs from `mqtl_susie_fine_mapping.sh`
7. Outputs from `pqtl_susie_fine_mapping.sh`

**Dependencies**: None

**Outputs**: None
