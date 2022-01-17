# Calderon *et al.* 2019

## Raw Data

The Calderon *et al.* study generated an atlas of transcriptomic and chromatin accessibility data for primary immune cell types in the blood under both unstimulated and stimulated conditions.

While 147 of the 166 RNA-Seq samples were run on an Illumina NovaSeq 6000 sequencer with 100 bp paired-end reads, 19 samples were run on an Illumina HiSeq 4000 sequencer with 76 bp paired-end reads.

Although there are 25 cell types present in the RNA-Seq samples, only 21 cell types have at least three samples (necessary for statistical tests).

Of the 8 donors for RNA-Seq samples, 4 made outsized contributions.

In the ATAC-Seq data, in contrast, 159 of the 175 samples were run on an Illumina HiSeq 4000 sequencer and 16 samples were run on an Illumina NovaSeq 6000 sequencer. All sequenced samples generated 76 bp paired-end reads.

Similar to the RNA-Seq data, 4 of the 7 sample donors made outsized contributions to the ATAC-Seq data.

| Covariate  | RNA-Seq | ATAC-Seq |
|:-----------|--------:|---------:|
| Donor      | 8       | 7        |
| Treatments | 2       | 2        |
| Lineages   | 6       | 6        |
| Cell Types | 25      | 25       |
| Sequencers | 2       | 2        |
| Reads      | Paired  | Paired   |
| Samples    | 166     | 175      |

**Note**: DCs are dendritic cells. pDCs are Plasmacytoid DCs.

| Lineage | Cell Type                | RNA-Seq Samples | ATAC-Seq Samples|
|:--------|:-------------------------|----------------:|----------------:|
| B       | Bulk_B                   | 7               | 7               |
| B       | Mem_B                    | 7               | 8               |
| B       | Naive_B                  | 8               | 7               |
| B       | Plasmoblasts             | 1               | 3               |
| CD8     | CD8pos_T                 | 8               | 7               |
| CD8     | Central_memory_CD8pos_T  | 7               | 8               |
| CD8     | Effector_memory_CD8pos_T | 8               | 8               |
| CD8     | Naive_CD8_T              | 8               | 8               |
| GD      | Gamma_delta_T            | 7               | 7               |
| CD4     | Effector_CD4pos_T        | 8               | 7               |
| CD4     | Follicular_T_Helper      | 8               | 9               |
| CD4     | Memory_Teffs             | 8               | 8               |
| CD4     | Memory_Tregs             | 7               | 8               |
| CD4     | Naive_Teffs              | 10              | 9               |
| CD4     | Regulatory_T             | 7               | 8               |
| CD4     | Th1_precursors           | 8               | 8               |
| CD4     | Th17_precursors          | 8               | 7               |
| CD4     | Th2_precursors           | 8               | 8               |
| CD4     | Naive_Tregs              | 5               | 4               |
| NK      | Immature_NK              | 1               | 5               |
| NK      | Mature_NK                | 10              | 10              |
| NK      | Memory_NK                | 1               | 6               |
| MYELOID | Monocytes                | 12              | 9               |
| MYELOID | Myeloid_DCs              | 3               | 3               |
| MYELOID | pDCs                     | 1               | 3               |

## Quality Control

### RNA-Seq

### ATAC-Seq

## Data Processing

### RNA-Seq

### ATAC-Seq

## Analysis

### RNA-Seq

### ATAC-Seq
