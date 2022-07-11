Analysis Data Folder
--------------------

To avoid exporting data to GitHub, the `data/` directory is ignored. Please add the following
directories to this `data/` folder. The scripts expect the following data files to also be
included within this folder.

## Expected Directory Structure

1. `./1000G/`
    1. `20131219.populations.tsv` - Retrieved from http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/ on November 2 2021
    2. `EUR.samples.txt` - Described below
    3. `integrated_call_samples_v3.20200731.ALL.ped` - Retrieved from http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ on November 2 2021 
2. `./Ensembl/`
    1. `GRCh38_to_GRCh37.chain.gz` - Retrieved from http://ftp.ensembl.org/pub/assembly_mapping/homo_sapiens/ on October 28 2021

Noting that the populations GBR, FIN, IBS, TSI, and CEU are in the EUR superpopulation, the following was used to generate the `EUR.samples.txt` file.

```
$ cat integrated_call_samples_v3.20200731.ALL.ped | grep "GBR\|FIN\|IBS\|TSI\|CEU" | awk '{ print $2; }' > EUR.samples.txt
```
