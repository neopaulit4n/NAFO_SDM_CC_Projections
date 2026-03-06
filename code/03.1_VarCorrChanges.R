
# Load variable correlation data ----
var_cor_df <- read_csv(list.files(file.path("output", "02_StepByStepOutputs"), pattern = "table_cor", full.names = TRUE),
  id = "iteration") %>%
  mutate(VME_group = gsub("^.+Outputs/(\\w+.+)_P\\d_.+$", "\\1", iteration),
         period = gsub("^.+_(P\\d)_.+$", "\\1", iteration),
         ssp = gsub("^.+_(\\d-\\d\\.\\d)_.+$", "\\1", iteration))
         

