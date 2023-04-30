# Mouse Mutant Cell Atlas (MMCA)

## Contents
- [Overview](#overview)
- [Scripts Download](#scripts-download)
- [Usage Instructions](#usage-instructions)
- [System Requirements](#system-requirements)
- [Demo](#demo)

# Overview

This repository contains data analysis scripts accompanying the manuscript "Single cell, whole embryo phenotyping of pleiotropic disorders of mammalian development", which is currently available on [bioRxiv] (https://www.biorxiv.org/content/10.1101/2022.08.03.500325v1). The comparitive analysis developed for MMCA can be found in section 4 (cell type composition analysis) and section 5 (lochNESS analysis). For more description on these analyses, please check out the methods section of the manuscript. 

# Scripts Download
The scripts and demo datasets can be downloaded through git: 
```
git clone https://github.com/shendurelab/MMCA.git
```
The download should only take a few seconds.

# Usage Instructions
These scripts are split into multiple sections corresponding to the major analysis in the order introduced in the manuscript. Running through the sccripts in the order provided will reproduce the major results in the manuscript. Unfortunately, the full dataset and some other intermediate files used in these scripts are too large to be hosted in this repo. Those data files and an interactive app to explore our dataset will be made freely available via https://atlas.gs.washington.edu/mmca_v2/. Please contact the authors if you find an intermedaite file missing.

# System requirements
The majority of these scripts were developed with `R/v4.0` on Linux or macOS operating systems. The dependencies of each script can be found at the beginning of each script. All dependencies can be dowloaded from `CRAN` or `BiocManager`.

# Demo
For those wanting to try it out for your own datasets, the best place to start is our demo scripts. For the lochNESS analysis, check out `demo_lochness.R`. The script takes a processed seurat object, and creates a dataframe with UMAP coordinates and lochNESS that can be used for plotting. In the demo, the dataset is seurat object containing Gli2 KO cells and wildtype cells in the Haematopoiesis trajectory.
