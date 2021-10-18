Analysis Data Folder
--------------------

To avoid exporting data to GitHub, the `data/` directory is ignored. Please add the following
directories to this `data/` folder. The scripts expect the following data files to also be
included within this folder.

## Expected Directory Structure

1. `./gene_coexpression_module_annotations/`
2. `./gene_expression/`
    1. `full-gains-SRS-predictions_mNN-RF.tsv`
    2. `Gene_info_864_20416.txt`
    3. `Logcpm_864_20416.txt`
    4. `Logcpm_864_20417_HLA.txt`
    5. `Sample_info_864.txt`
    6. `Sample_key.txt`
3. `./gene_expression_generated/`
4. `./protein_expression/`
    1. `data_291x1860_MS2019.csv`
    2. `protein_info_291_MS2019.csv`
    3. `sample_info_1860_MS2019.csv`
5. `./protein_expression_generated/`
