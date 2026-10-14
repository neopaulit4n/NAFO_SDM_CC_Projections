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

- **R**: all the code, in three subfolders. Function definitions, once
  extracted from the scripts, go directly in R/ (top level).
  - **analysis**: the scripts that produce the results.
    - 00_CodeCaller.R: Main script that runs the other scripts in this
      folder in turn and writes all model outputs to the **output**
      folder (below). Run it from the project root
      (`source("R/analysis/00_CodeCaller.R")`), after the scripts in
      **R/data-raw**.
    - 01\_ to 12\_: called (“sourced”) from and described in
      00_CodeCaller.R.
    - 13\_ to 17\_: not called by 00_CodeCaller.R; run them afterwards
      in the same R session (they use its objects), for specific VME
      groups or figures/tables.
  - **data-raw**: scripts used to clean the original input data. They
    only need to be run once; the cleaned data are written to
    data/processed and then used by the modelling scripts.
  - **archive**: past exploratory and reference code, kept for the record
    and not used by the pipeline (see R/archive/README.md).
- **data**
  - **raw**: raw, unedited data
  - **processed**: cleaned, edited data obtained using the R/data-raw
    scripts from the raw data
- **output**: model outputs automatically get stored here, stored in
  separate subfolders for each VME group.
- **doc**: reports, including the reproducibility review
  (doc/reproducibility.qmd).
