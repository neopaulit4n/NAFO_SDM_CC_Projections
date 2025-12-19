
# Editing and transferring SDM2024 code for use with CC projection model data

library(tidyverse)

# Load data ----
resp_df <- read_csv("data/raw/SDM2024/FINAL Response Variables/Black corals/black_corals.csv",
                    show_col_types = FALSE)

