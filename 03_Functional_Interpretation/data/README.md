Analysis Data Folder
--------------------

To avoid exporting data to GitHub, the `data/` directory is ignored. Please add the following
directories to this `data/` folder. The scripts expect the following data files to also be
included within this folder.

## Expected Directory Structure

1. `hg38-blacklist.v2.bed` - Retrieved from [https://github.com/Boyle-Lab/Blacklist](https://github.com/Boyle-Lab/Blacklist)
2. `41588_2019_505_MOESM6_ESM` - Retrieved from [https://www.nature.com/articles/s41588-019-0505-9#Sec30](https://www.nature.com/articles/s41588-019-0505-9#Sec30) - Supplementary Table 3
3. The `collapse_annotation.py` script from GTEx (Retrieved from [https://github.com/broadinstitute/gtex-pipeline/tree/master/gene_model](https://github.com/broadinstitute/gtex-pipeline/tree/master/gene_model)) was used to collapse the Ensembl 99 annotation file. This generated the `Homo_sapiens.GRCh38.99.collapsed.gtf` file.
