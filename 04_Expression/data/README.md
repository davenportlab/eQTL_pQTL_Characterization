Analysis Data Folder
--------------------

To avoid exporting data to GitHub, the `data/` directory is ignored. Please add the following
directories to this `data/` folder. The scripts expect the following data files to also be
included within this folder.

## Expected Directory Structure

1. `./gene_expression/`
    1. `Cell_props_864.txt` - Cell type proportions for RNA-seq samples
2. `xCell_Aran_et_al_Additional_File_1.csv` - Retrieved from [https://genomebiology.biomedcentral.com/articles/10.1186/s13059-017-1349-1](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-017-1349-1). Converted to CSV format using Excel after correction of minor spelling errors.
3. `gwas_catalog_v1.0.2-studies_r2022-02-21.tsv` - Retrieved from [https://www.ebi.ac.uk/gwas/docs/file-downloads](https://www.ebi.ac.uk/gwas/docs/file-downloads). Downloaded on 23 February 2022.
4. `gwas_catalog_v1.0.2-associations_e105_r2022-02-21.tsv` - Retrieved from [https://www.ebi.ac.uk/gwas/docs/file-downloads](https://www.ebi.ac.uk/gwas/docs/file-downloads). Downloaded on 23 February 2022.
5. `clindata_2022-06-21_11_08-38.csv` - Clinical Data retrieved from Sepsis LabKey server
6. `out_2022-06-21_11-16-02.csv` - Outcome Data retrieved from Sepsis LabKey server
7. `kwok_2022_cell_markers.csv` - Data retrieved from supplementary files of manuscript from [https://www.medrxiv.org/content/10.1101/2022.03.22.22272723v1](https://www.medrxiv.org/content/10.1101/2022.03.22.22272723v1)
8. `kwok_2022_cell_types.csv` - Data retrieved from supplementary files of manscript from [https://www.medrxiv.org/content/10.1101/2022.03.22.22272723v1](https://www.medrxiv.org/content/10.1101/2022.03.22.22272723v1)
