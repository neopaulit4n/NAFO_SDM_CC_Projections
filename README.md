# Code repository for the NAFO Regulatory Area Species Distribution Modelling under Climate Change Projection Scenarios


# Summary

This R code repository contains the code necessary to run the NAFO
Regulatory Area species distribution models under climate change
projection scenarios.

The following NAFO SCRs/publications supported by these models are: TBD

# Input data

Data used for running these models can be found here: TBD

# Usage guide and structure

The project structure is as follows:

- **code**
  - **01_DataCleanup**: Subfolder which contains scripts used to clean
    the original input data. These scripts only need to be ran once and
    the intermediary cleaned input data get automatically placed in the
    data/processed folder and are then used for the modelling scripts.
  - **02_Exploratory**: Past exploratory code used to test different
    things before a set protocol was established; only used as reference
    code.
  - 00_CodeCaller.R: Main script that runs all the other scripts in this
    folder in turn to output all relevant model outputs into the
    **output** folder (below). You will only need to run this script to
    output everything after running the scripts from the 01_DataCleanup
    folder.
  - …. All other scripts in this folder are called (“sourced”) from and
    described in the 00_CodeCaller.R script.
- **data**
  - **raw**: raw, unedited data
  - **processed**: cleaned, edited data obtained using the
    01_DataCleanup scripts from the raw data
- **output**: model outputs automatically get stored here, stored in
  separate subfolders for each VME group.
