# pQTL Mapping

This folder contains the scripts for processing the pQTL mapping results. The order of running the scripts is described below:

1. [pQTL Processing](#pqtl-processing)
    1. [Protein Metadata](#protein-metadata)
    2. [Extract pQTL](#extract-pqtl)
    3. [pQTL Analysis](#pqtl-analysis)

---

## pQTL Processing

Run the scripts in the following order:

---

### Protein Metadata

**File**: `protein_metadata.ipynb`

**Description**: Generates metadata for the proteins assayed using mass spectrometry.

**Expected Inputs**:
1. UniProt Data
2. Protein Information from Mass Spectrometry
3. Gene Information from RNA-seq

**Dependencies**: None

**Outputs**:
1. Table of Metadata (Gene-to-Protein Mapping)

---

### Extract pQTL

**File**: `extract_pqtl.sh`

**Description**: Take raw outputs from QTL mapping performed at WHI and generate clean data.

**Expected Inputs**:
1. Outputs from `protein_metadata.ipynb`
2. Patient Genotyped Variants
3. pQTL Summary Statistics

**Dependencies**:
1. `extract_pqtl.nf`
2. `extract_pqtl.R`
3. `extract_pqtl_aggregate.R`

**Outputs**:
1. cis-pQTL
2. trans-pQTL

---

### pQTL Analysis

**File**: `pqtl.ipynb`

**Description**: Identify significant cis-pQTL and trans-pQTL.

**Expected Inputs**:
1. Outputs from `protein_metadata.ipynb`
2. Outputs from `extract_pqtl.sh`

**Dependencies**: None

**Outputs**:
1. Significant cis-pQTL
2. Significant trans-pQTL
