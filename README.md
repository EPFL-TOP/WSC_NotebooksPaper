# Repository Overview

This repository contains annotated Python notebooks, along with the in vitro and in vivo data used in the workflow to produce the figures and results in our eLife paper:
**“Cell-autonomous timing drives the vertebrate segmentation clock’s wave pattern”**
 
[DOI: https://doi.org/10.7554/eLife.93764](https://doi.org/10.7554/eLife.93764).

Our raw imaging data can be downloaded at https://sv-open.epfl.ch/upoates-public/Rohde_Bercowsky_eLife2024_Data/

Given the nature of our workflow, the data, analysis code, and results contributing to a given eLife figure are distributed throughout Jupyter notebooks organized by subject. Please start by navigating to the folder **“notebooks”**. The exact data used in each analysis is specified in the respective notebook and can be found in the **“Data”** folder. 

> **Note:** Any figure numbers in the Data folder do not correspond directly to the paper.

To ensure the notebooks run without issues, please maintain the repository's directory structure after download. Results in the form of PDF plots are saved inside the **“Results”** folder located within the **“Data”** folder, and have descriptive titles.

## Environment Setup

To run the notebooks, please install the following conda environment:

```bash
conda env create -f wsc.yml
```

# Table of Contents of Notebooks

DOI:[https://doi.org/10.7554/eLife.93764](https://doi.org/10.7554/eLife.93764)

The analysis pipeline connecting experimental data and results for each figure in Rohde et al., 2024 can be found in the following Jupyter notebooks:

### eLife Figure 1 and Figure 1 supplements:
- PSM4 Dynamics in the embryo Part I
- PSM4 Dynamics in the embryo Part II
- PSM4 Dynamics in culture Part I
- PSM4 Dynamics in culture Part II
- Dynamics by cycle in culture
- Dynamics by cycle in the embryo

### eLife Figure 2 and Figure 2 supplements:
- Backtracked somites
- PSM4 Dynamics in culture Part II
- PSM4 Dynamics in the embryo Part I
- PSM4 Dynamics in the embryo Part II

### eLife Figure 3 and Figure 3 supplements:
- FGF in culture Mesp analysis
- FGF in culture Her1 analysis

### eLife Figure 4 and Figure 4 supplements:
- TB dynamics in embryo
- AP Dynamics in culture 1 & 2
