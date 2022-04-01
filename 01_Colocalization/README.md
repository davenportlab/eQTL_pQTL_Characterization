# Colocalization Workflow

1. `scripts/gains_eqtl_data_preparation.R` - Prepare eQTL data for conversion from GRCh38 to GRCh37 
coordinates. Save standard deviation of expression data.
2. `workflow/liftover_snp_coordinates.sh` - Use chain file to convert coordinates between genome
builds.
3. `scripts/gains_eqtl_data_preparation_hg19_coordinates.R` - Split the cis-eQTL by the cognate 
gene. Attach the GRCh37 coordinates to each locus.
4. `scripts/gains_pilot_pqtl_data_preparation.R` - Split the cis-pQTL by the cognate protein.
Calcualte standard deviation of the expression data.
5. `scripts/coloc_ocv_assumption.R` - Run the original version of COLOC with sensitivity analysis.
6. `scripts/coloc_susie.R` - Run COLOC with SuSiE.

## Extraction Plan

1. R to convert rsIDs to position/base format
2. R to prepare eQTL and pQTL data (no need to convert to GRCh37)
3. Python to run VCF tools to extract LD values online
4. COLOC with LD Matrix