Analysis Data Folder
--------------------

To avoid exporting data to GitHub, the `data/` directory is ignored. Please add the following
directories to this `data/` folder. The scripts expect the following data files to also be
included within this folder.

## Expected Directory Structure

1. `hg38-blacklist.v2.bed` - Retrieved from [https://github.com/Boyle-Lab/Blacklist](https://github.com/Boyle-Lab/Blacklist)
2. `41588_2019_505_MOESM6_ESM` - Retrieved from [https://www.nature.com/articles/s41588-019-0505-9#Sec30](https://www.nature.com/articles/s41588-019-0505-9#Sec30) - Supplementary Table 3
3. `41588_2019_505_MOESM8_ESM` - Retrieved from [https://www.nature.com/articles/s41588-019-0505-9#Sec30](https://www.nature.com/articles/s41588-019-0505-9#Sec30) - Supplementary Table 5
4. `Ram_Mohan_et_al_2022_Table_S1.xlsx` - Retrieved from [https://www.life-science-alliance.org/content/4/8/e202000976/tab-figures-data](https://www.life-science-alliance.org/content/4/8/e202000976/tab-figures-data) - Supplementary Table 1
5. The `collapse_annotation.py` script from GTEx (Retrieved from [https://github.com/broadinstitute/gtex-pipeline/tree/master/gene_model](https://github.com/broadinstitute/gtex-pipeline/tree/master/gene_model)) was used to collapse the Ensembl 99 annotation file. This generated the `Homo_sapiens.GRCh38.99.collapsed.gtf` file.
6. `hg19ToHg38.over.chain` - Retrieved from [http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/](http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/). Used `gunzip -c hg19ToHg38.over.chain.gz` to unzip the file.
7. `GSE118189_ATAC_counts.txt.gz` - Retrieved from [https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE118189](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE118189).
8. `JASPAR2022_CORE_vertebrates_non-redundant_pfms_meme.txt` - Retrieved from [https://jaspar.genereg.net/download/data/2022/CORE/JASPAR2022_CORE_vertebrates_non-redundant_pfms_meme.txt](https://jaspar.genereg.net/download/data/2022/CORE/JASPAR2022_CORE_vertebrates_non-redundant_pfms_meme.txt).