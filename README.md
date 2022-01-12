# eQTL pQTL Characterization

This is part of Nikhil's Thesis Project for the 2021/2022
academic year.

## Summary

The goal of this project is to characterize mechanisms
underlying molecular QTL in the Genomic Advances in 
Sepsis ([GAinS](https://ukccggains.com/)) cohort.

## Technical Details

We use [Anavonda environments](https://docs.conda.io/projects/conda/en/latest/index.html)
to manage dependencies. Most analyses are run using
either R or Python. Analytical steps were coordinated on 
a computing cluster using
[Nextflow](https://www.nextflow.io/).

A minimal setup for running the following analyses is a
Linux workstation or compute cluster with:

1. An installation of Anaconda
2. An [executor](https://www.nextflow.io/docs/latest/executor.html)
    supported by Nextflow

## Analyses

The project is divided into five broad analytical steps.
Each analysis has a dedicated root-level directory. Each
step also has its own Anaconda environment. A project-level
environment is used for the Jupyter Notebooks.

### Analysis Steps

1. `01_Colocalization/` - The goal is to colocalize the
    eQTL from the GAinS cohort with pQTL from GAinS and
    other sepsis-related cohorts.
2. `02_pQTL_Mapping` - The goal is to process pQTL data
    generated by the Knight group at Oxford.
3. `03_Functional_Interpretation` - The goal is to overlay
    datasets to better understand the mechanisms
    underlying eQTL detected in GAinS.
4. `04_Expression` - The goal is to jointly analyze gene
    and protein expression.
5. `05_MR` - The goal is to perform Mendelian
    Randomization to hypothesize about causality between
    molecular exposures and outcomes after the
    colocalization analysis.

## Setup

### Locate Required Data

For some steps, external data is required. This data is
not uploaded to the repository. Details on how to acquire
this data is included in a `README.md` file in the `data/`
directory for each step.

### Install Anaconda Environments

Begin by installing the Anaconda environments for the
project. The project-level Jupyter environment is stored
in `env/` and the environment for each step is stored in
the `env/` directory of the step.

As an example, to install the Jupyter environment, start by
checking if Anaconda is installed.

```
$ conda --version
conda 4.9.2
```

Install the environment to some location.

```
$ conda env create \
>   --file jupyter_env/environment.yml \
>   --prefix path/to/install/jupyter_env/
```

Activate the environment when required.

```
$ conda activate path/to/install/jupyter_env/
(jupyter_env) $ 
```

### Detailed Steps for Replication

The steps required to replicate analyses are included in each step's directory. Please refer to the `README.md` file
included for each step.
