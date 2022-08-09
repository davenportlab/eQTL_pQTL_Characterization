# Data Management

## Sample Metadata

Sample SRA accessions for the ATAC-seq and RNA-seq studies are stored in this folder.

1. `reads_atac_seq.txt` - Contains sample information for each ATAC-seq run
2. `reads_rna_seq.txt` - Contains sample information for each RNA-seq run
3. `eQTL_pQTL_metadata.tsv` - Contains mapping between eGenes and pGenes

## Terminal Data Storage

> Keep on Lustre

**Directory**: `/lustre/scratch119/humgen/projects/gains_team282/nikhil/colocalization/cis_eqtl/conditional_effects/LMM`

**Rationale**: These are conditional summary statistics derived from the conditional cis-eQTL analysis. In brief, for each individual signal, lead eSNPs from all other conditional cis-eQTL for the eGene were added to the LMM.

> Keep on Lustre

**Directory**: `/lustre/scratch119/humgen/projects/gains_team282/nikhil/colocalization/cis_eqtl/fine_mapping`

**Rationale**: Contains all the fine mapping results for cis-eQTL. This includes FINEMAP and SuSiE run on the original summary statistics and the conditional summary statistics. It also includes LD "tagging SNP sets" for all the molecular QTL.

> Keep on Lustre

**Directory**: `/lustre/scratch119/humgen/projects/gains_team282/nikhil/colocalization/mqtl/fine_mapping`

**Rationale**: Contains all the fine mapping results for the module QTL (FINEMAP and SuSiE).

> Keep on Lustre

**Directory**: `/lustre/scratch119/humgen/projects/gains_team282/nikhil/colocalization/pqtl/fine_mapping`

**Rationale**: Contains all the fine mapping results for the pQTL (FINEMAP and SuSiE).

> Keep on Lustre

**Directory**: `/lustre/scratch119/humgen/projects/gains_team282/nikhil/data/GAinS_Microarray`

**Rationale**: This directory contains the `.rds` object provided by Eddie. It contains the microarray expression data from GAinS samples, collapsed by genes.

> Keep on Lustre

**Directory**: `/lustre/scratch119/humgen/projects/gains_team282/nikhil/data/SNPsnap`

**Rationale**: This directory contains results from the SNPsnap web server. This includes 10,000 SNPs matched to the lead conditional cis-eQTL, which can be useful for enrichment testing.

> Keep on Lustre

**Directory**: `/lustre/scratch119/humgen/projects/gains_team282/epigenetics/accessibility/merged`

**Rationale**: This directory contains merged alignments from the ATAC-seq and RNA-seq data. This includes the BAM files after alignment, peaks called by MACS2 for the ATAC-seq data, bigWig and bedGraph files for viewing read pileups, and RNA-seq gene counts.

> Archive

**Directory**: `/lustre/scratch119/humgen/projects/gains_team282/nikhil/expression`

**Rationale**: This directory contains all the results from the WGCNA analysis. This includes module generation, module annotation, and module QTL steps of the analysis. Individual files mentioned in the data management spreadsheet can be retained on Lustre, but in general it will be beneficial to archive the rest.

> Archive

**Directory**: `/lustre/scratch119/humgen/projects/gains_team282/epigenetics/vep`

**Rationale**: This directory contains results from the VEP analysis on the lead conditional cis-eQTL. VCF files were generated using [this script](https://github.com/davenportlab/eQTL_pQTL_Characterization/tree/main/03_Functional_Interpretation#preparation) and analysed using [this script](https://github.com/davenportlab/eQTL_pQTL_Characterization/tree/main/03_Functional_Interpretation#ensembls-vep) and [this script](https://github.com/davenportlab/eQTL_pQTL_Characterization/tree/main/03_Functional_Interpretation#vep-analysis). Results were not exactly useful for interpretation, but may be good for querying individual variants.

> Archive

**Directory**: `/lustre/scratch119/humgen/projects/gains_team282/epigenetics/enrichment/cheers`

**Rationale**: Contains results from running CHEERS in the neutrophil data. The results are challenging to interpret without more conditions and better quality data. The data is generated using [`cheers.sh`](https://github.com/davenportlab/eQTL_pQTL_Characterization/blob/main/03_Functional_Interpretation/workflow/enrichment/cheers.sh), detailed [here](https://github.com/davenportlab/eQTL_pQTL_Characterization/tree/main/03_Functional_Interpretation#cheers)

> Archive

**Directory**: `/lustre/scratch119/humgen/projects/gains_team282/epigenetics/enrichment/go_shifter`

**Rationale**: Contains results from running GoShifter on the immune and neutrophil datasets. The data is generated using [`go_shifter.sh`](https://github.com/davenportlab/eQTL_pQTL_Characterization/blob/main/03_Functional_Interpretation/workflow/enrichment/go_shifter.sh), detailed [here](https://github.com/davenportlab/eQTL_pQTL_Characterization/tree/main/03_Functional_Interpretation#goshifter).

> Archive

**Directory**: `/lustre/scratch119/humgen/projects/gains_team282/epigenetics/enrichment/partitioned_heritability`

**Rationale**: A partitioned heritability model was used to estimate module QTL heritability in various genomic annotations. The data is generated using [`partitioned_heritability.sh`](https://github.com/davenportlab/eQTL_pQTL_Characterization/blob/main/03_Functional_Interpretation/workflow/enrichment/partitioned_heritability.sh), detailed [here](https://github.com/davenportlab/eQTL_pQTL_Characterization/tree/main/03_Functional_Interpretation#partitioned-heritability).

> Archive

**Directory**: `/lustre/scratch119/humgen/projects/gains_team282/epigenetics/enrichment/snp_snap`

**Rationale**: Here, I tested matched SNPs and the cis-eQTL for overlap with various genomic annotations. This was used to generate the matched SNP enrichment results in my thesis. The data is generated using [`snp_snap_matching.sh`](https://github.com/davenportlab/eQTL_pQTL_Characterization/blob/main/03_Functional_Interpretation/workflow/enrichment/snp_snap_matching.sh), detailed [here](https://github.com/davenportlab/eQTL_pQTL_Characterization/tree/main/03_Functional_Interpretation#matched-snp-enrichment).

> Archive

**Directory**: `/lustre/scratch119/humgen/projects/gains_team282/epigenetics/accessibility/processed/atac_seq/multiqc`

**Rationale**: Quality control metrics from the alignment of the publicly available ATAC-seq data.

> Archive

**Directory**: `/lustre/scratch119/humgen/projects/gains_team282/epigenetics/accessibility/processed/rna_seq/multiqc`

**Rationale**: Quality control metrics from the alignment of the publicly available RNA-seq data.

> Archive

**Directory**: `/lustre/scratch119/humgen/projects/gains_team282/epigenetics/accessibility/analysis`

**Rationale**: This directory contains all the analyses conducted using the RNA-seq and ATAC-seq data. This includes peak sets, TSS enrichment scores, and motif enrichment analysis. It might be useful to keep certain outputs on Lustre that are mentioned in the data management spreadsheet, but otherwise the entire directory can be archived.

> Delete

**Directory**: `/nfs/users/nfs_n/nm18/eQTL_pQTL_Characterization`

**Reproduce**: This is the GitHub repository. It can be cloned from [here](https://github.com/davenportlab/eQTL_pQTL_Characterization). Instructions for finding and adding data to the repository directory structure is included in `README.md` files in each main folder.

> Delete

**Directory**: `/lustre/scratch119/humgen/projects/gains_team282/nikhil/bin`

**Reproduce**: This directory contains binaries and scripts that could not be installed into the `conda` environment. This includes the MEME suite, FINEMAP, RNASeQC, [CHEERS](https://github.com/TrynkaLab/CHEERS/tree/python3), and my version of [GoShifter](https://github.com/NMilind/LeanGoShifter). The first three are discussed in the main `README.md` file. The last two can be cloned from their respective repositories.

> Delete

**Directory**: `/lustre/scratch119/humgen/projects/gains_team282/nikhil/colocalization/cis_eqtl/conditional_effects/COJO`

**Reproduce**: This is an old, experimental analysis that was performed using GCTA's [COJO](https://yanglab.westlake.edu.cn/software/gcta/#COJO). The idea was to to derive conditional summary statistics for each signal. However, I have gone back and run our original LMM to derive conditional summary statistics using the original genotyping data. There is no script on the repository that can regenerate this data.

> Delete

**Directory**: `/lustre/scratch119/humgen/projects/gains_team282/nikhil/colocalization/cis_eqtl/conditional_effects/FINEMAP`

**Reproduce**: This is an old, experimental analysis that was performed by first using FINEMAP on the summary statistics and using the lead eSNP in each credible set to derive conditional summary statistics. However, I have gone back and run our original LMM to derive conditional summary statistics using the original genotyping data. There is no script on the repository that can regenerate this data.

> Delete

**Directory**: `/lustre/scratch119/humgen/projects/gains_team282/nikhil/colocalization/pqtl/conditional_effects`

**Reproduce**: This is an old, experimental analysis that was performed by first using FINEMAP on the summary statistics and using the lead eSNP in each credible set to derive conditional summary statistics. There is no script on the repository that can regenerate this data.

> Delete

**Directory**: `/lustre/scratch119/humgen/projects/gains_team282/nikhil/conda`

**Reproduce**: This is the `conda` environment. The repository contains two methods to reproduce this environment. The [`environment.yaml`](https://github.com/davenportlab/eQTL_pQTL_Characterization/blob/main/env/environment.yml) file can be used to create the `conda` environment based on a subset of required software. The asset in the [Thesis Submission release](https://github.com/davenportlab/eQTL_pQTL_Characterization/releases/tag/Thesis-Submission) contains a dump of all the software in this `conda` environment, which can also be used to create a new environment.

> Delete

**Directory**: `/lustre/scratch119/humgen/projects/gains_team282/nikhil/miniconda3`

**Reproduce**: I installed [`miniconda3`](https://docs.conda.io/en/latest/miniconda.html) to run `conda`. Eventually, I started using [`mamba`](https://github.com/mamba-org/mamba) instead.

> Delete

**Directory**: `/lustre/scratch119/humgen/projects/gains_team282/nikhil/data/1000G`

**Reproduce**: This directory contains genotypes from the 1000 Genomes project at module QTL loci. This was going to be used to perform fine mapping on the GWAS association statistics before performing colocalisation with the module QTL. This data can be regenerated using the [`extract_mqtl_ld_1000g.sh`](https://github.com/davenportlab/eQTL_pQTL_Characterization/blob/main/01_Colocalization/workflow/SuSiE_fine_mapping/extract_mqtl_ld_1000g.sh) script, detailed [here](https://github.com/davenportlab/eQTL_pQTL_Characterization/tree/main/01_Colocalization#extract-1000g-ld).

> Delete

**Directory**: `/lustre/scratch119/humgen/projects/gains_team282/nikhil/data/COMBAT`

**Reproduce**: Contains publicly-available data from the COMBAT paper. The data was downloaded from [here](https://zenodo.org/record/6120249).

> Delete

**Directory**: `/lustre/scratch119/humgen/projects/gains_team282/nikhil/data/EBI_GWAS_Catalog`

**Reproduce**: The accessions of the GWAS studies used for colocalisation is in my thesis. Summary statistics can be downloaded from the [EBI GWAS Catalog](https://www.ebi.ac.uk/gwas/downloads/summary-statistics). I also converted any coordinates from hg19 to GRCh38.

> Delete

**Directory**: `/lustre/scratch119/humgen/projects/gains_team282/nikhil/data/ENCODE`

**Reproduce**: I downloaded all human cCREs from the ENCODE project [website](https://screen.encodeproject.org/) under the Downloads tab. The downloaded file was GRCh38-cCREs.bed. I used the following command to split the BED file based on the cCRE type.

```
awk \
        'BEGIN { OFS="\t"; }
        { 
                if ($6 ~ "CTCF-bound") {
                        x = $6; 
                        gsub(",.*", "", x);
                        print $1, $2, $3, $4, $5, x, "\tCTCF-bound" > "GRCh38-cCREs." x ".bed";
                } else {
                        print $1, $2, $3, $4, $5, x, "\tCTCF not bound" > "GRCh38-cCREs." x ".bed";
                }
        }' GRCh38-cCREs.bed
```

> Delete

**Directory**: `/lustre/scratch119/humgen/projects/gains_team282/nikhil/data/Expecto`

**Reproduce**: This is an old, experimental analysis. I ran all lead molecular QTL through [Expecto](https://hb.flatironinstitute.org/expecto/).

> Delete

**Directory**: `/lustre/scratch119/humgen/projects/gains_team282/nikhil/data/genotypes`

**Reproduce**: This directory contains any genotypes that were extracted from the binary PLINK format files at specific loci. These were generally used to perform single-variant association mapping and to generate LD matrices from the cohort. They can be regenerated using multiple scripts:

1. [`extract_genotypes.sh`](https://github.com/davenportlab/eQTL_pQTL_Characterization/blob/main/01_Colocalization/workflow/extract_genotypes.sh), detailed [here](https://github.com/davenportlab/eQTL_pQTL_Characterization/tree/main/01_Colocalization#extract-genotypes)
2. [`extract_genotypes.sh`](https://github.com/davenportlab/eQTL_pQTL_Characterization/blob/main/04_Expression/workflow/extract_genotypes.sh), detailed [here](https://github.com/davenportlab/eQTL_pQTL_Characterization/tree/main/04_Expression#prepare-genotypes-for-module-qtl)

> Delete

**Directory**: `/lustre/scratch119/humgen/projects/gains_team282/nikhil/data/ROADMAP`

**Reproduce**: This directory contains select Roadmap epigenomes, accessions for which are catalogued in my thesis. Briefly, using the Roadmap epigenomes [metadata](https://github.com/davenportlab/eQTL_pQTL_Characterization/blob/main/03_Functional_Interpretation/metadata/roadmap_epigenomics_chromhmm.csv), I ran the following script:

```
meta=/nfs/users/nfs_n/nm18/eQTL_pQTL_Characterization/03_Functional_Interpretation/metadata/roadmap_epigenomics_chromhmm.csv

awk -F ',' 'NR  1 { print $4; }' $meta | xargs wget

ls -1 | grep gz | xargs gunzip

for bed_file in E*_18_core_K27ac_hg38lift_segments.bed
do

    awk '{
        for (i = 1; i <= 18; i++) {
            regex = "E"i;
            if ($4 ~ regex) {
                gsub("\\.bed", "", FILENAME);
                print $0 > FILENAME "-E" i ".bed";
            }
        }
    }' $bed_file

done
```

> Delete

**Directory**: `/lustre/scratch119/humgen/projects/gains_team282/nikhil/data/VEP`

**Reproduce**: This was an install of VEP that I was attempting. In the end, I used Sanger Institute's in-house VEP install. This directory is thus obsolete.

> Delete

**Directory**: `/lustre/scratch119/humgen/projects/gains_team282/nikhil/data/Yazar_et_al_2022`

**Reproduce**: Processed scRNA-seq results from [Yazar *et al.*](https://www.science.org/doi/10.1126/science.abf3041) in both Seurat and H5AD format was retrieved from [here](https://cellxgene.cziscience.com/collections/436154da-bcf1-4130-9c8b-120ff9a888f2).

> Delete

**Directory**: `/lustre/scratch119/humgen/projects/gains_team282/nikhil/functional_interpretation`

**Reproduce**: I overlapped lead and conditional eSNPs with DA peaks and caQTL reported by [Calderon *et al.*](https://www.nature.com/articles/s41588-019-0505-9) in their supplementary data. There is no script in the repository that reproduces this analysis, although [this script](https://github.com/davenportlab/eQTL_pQTL_Characterization/blob/main/03_Functional_Interpretation/scripts/atac_08_calderon_et_al_comparison.ipynb) can be used as a starting point to work with the supplementary data.

> Delete

**Directory**: `/lustre/scratch119/humgen/projects/gains_team282/epigenetics/bowtie2_genome_index`

**Reproduce**: This is a Bowtie 2 index for read alignment. It was based on HGI's version of GRCh38 (Ensembl version 99). It can be generated using [`index_genome.sh`](https://github.com/davenportlab/eQTL_pQTL_Characterization/blob/main/03_Functional_Interpretation/workflow/index_genome.sh), detailed [here](https://github.com/davenportlab/eQTL_pQTL_Characterization/tree/main/03_Functional_Interpretation#index-genomes).

> Delete

**Directory**: `/lustre/scratch119/humgen/projects/gains_team282/epigenetics/star_genome_index`

**Reproduce**: This is a STAR index for read alignment. It was based on HGI's version of GRCh38 (Ensembl version 99). It can be generated using [`index_genome.sh`](https://github.com/davenportlab/eQTL_pQTL_Characterization/blob/main/03_Functional_Interpretation/workflow/index_genome.sh), detailed [here](https://github.com/davenportlab/eQTL_pQTL_Characterization/tree/main/03_Functional_Interpretation#index-genomes).

> Delete

**Directory**: `/lustre/scratch119/humgen/projects/gains_team282/epigenetics/calderon_et_al_hg19`

**Reproduce**: This is an old, experimental analysis where I was attempting to use CHEERS and GoShifter on the supplementary data provided by [Calderon *et al.*](https://www.nature.com/articles/s41588-019-0505-9). The CHEERS analysis can be reproduced using [`cheers_calderon_et_al_hg19`](https://github.com/davenportlab/eQTL_pQTL_Characterization/blob/main/03_Functional_Interpretation/workflow/calderon_et_al_hg19/cheers_calderon_et_al_hg19.sh), detailed [here](https://github.com/davenportlab/eQTL_pQTL_Characterization/tree/main/03_Functional_Interpretation#cheers-on-calderon-et-al-data). The  GoShifter analysis can be reproduced using [`goshifter_calderon_et_al_hg19.sh`](https://github.com/davenportlab/eQTL_pQTL_Characterization/blob/main/03_Functional_Interpretation/workflow/calderon_et_al_hg19/goshifter_calderon_et_al_hg19.sh), detailed [here](https://github.com/davenportlab/eQTL_pQTL_Characterization/tree/main/03_Functional_Interpretation#goshifter-on-calderon-et-al-data).

> Delete

**Directory**: `/lustre/scratch119/humgen/projects/gains_team282/epigenetics/regulation`

**Reproduce**: This is an experimental analysis where I was trying to identify peaks that changed shape upon stimulation. This can be reproduced using [`peak_shape_features.sh`](https://github.com/davenportlab/eQTL_pQTL_Characterization/blob/main/03_Functional_Interpretation/workflow/regulation/peak_shape_features.sh), detailed [here](https://github.com/davenportlab/eQTL_pQTL_Characterization/tree/main/03_Functional_Interpretation#peak-shapes).

> Delete

**Directory**: `/lustre/scratch119/humgen/projects/gains_team282/epigenetics/enrichment/heritability`

**Reproduce**: The partitioned heritability model required genetic relationship matrices (GRMs), which were generated for various annotations using GCTA-GREML. This can be reproduced using [`heritability.sh`](https://github.com/davenportlab/eQTL_pQTL_Characterization/blob/main/03_Functional_Interpretation/workflow/enrichment/heritability.sh), detailed [here](https://github.com/davenportlab/eQTL_pQTL_Characterization/tree/main/03_Functional_Interpretation#create-grms).

> Delete

**Directory**: `/lustre/scratch119/humgen/projects/gains_team282/epigenetics/enrichment/snp_matching`

**Reproduce**: Before using SNPsnap, I attempted to generate sets of matched SNPs using the cohort's genotyping data. This data can be reproduced using [`snp_matching.sh`](https://github.com/davenportlab/eQTL_pQTL_Characterization/blob/main/03_Functional_Interpretation/workflow/enrichment/snp_matching.sh), detailed [here](https://github.com/davenportlab/eQTL_pQTL_Characterization/tree/main/03_Functional_Interpretation#snp-matching-within-cohort).

> Delete

**Directory**: `/lustre/scratch119/humgen/projects/gains_team282/epigenetics/accessibility/raw`

**Reproduce**: This directory contains raw FASTQ files from publicly-available RNA-seq and ATAC-seq experiments. The data can be reproduced using [`download_from_sra.sh`](https://github.com/davenportlab/eQTL_pQTL_Characterization/blob/main/03_Functional_Interpretation/workflow/download_from_sra.sh), detailed [here](https://github.com/davenportlab/eQTL_pQTL_Characterization/tree/main/03_Functional_Interpretation#download-from-sra).

> Delete

**Directory**: `/lustre/scratch119/humgen/projects/gains_team282/epigenetics/accessibility/processed`

**Reproduce**: This directory contains trimmed FASTQ files and alignments for the ATAC-seq and RNA-seq data. The data generation is a bit complex, and is detailed fully [here](https://github.com/davenportlab/eQTL_pQTL_Characterization/tree/main/03_Functional_Interpretation).
