
# Reading in table of variable importance created manually from 2024/2025 NAFO SCR BNAM SDM papers
# Mean decrease in GINI is approximated from the graphs included in those papers since the original numbers were not available

var_imp <- read.csv("data/processed/Approx_2024_2025_SCR_VarImportance.csv") %>%
  # Calculate difference in mean GINI decrease between each row by VME group adjusted by max mean decrease in GINI for that group
  group_by(VME_Group) %>%
  mutate(max_gini = max(mean_dec_GINI)) %>%
  ungroup()

# Take the first or if tied first two top static variables based on mean decrease in GINI for each functional group/taxon
static_var_top <- var_imp %>%
  filter(var_type == "static") %>%
  group_by(VME_Group) %>%
  mutate(gini_diff = (mean_dec_GINI - lag(mean_dec_GINI))/max_gini*-100,
         lagged_gini_diff = lead(gini_diff))

static_var_top_f <- static_var_top %>%
  filter(row_number() == 1 | (row_number() == 2 & gini_diff <= 5)) %>%
  ungroup() %>%
  select(VME_Group, variable, mean_dec_GINI, lagged_gini_diff)

write.csv(static_var_top_f, 
          "data/processed/TopStaticVars_ByVMEGroup_2024_2025_SCR.csv", row.names = FALSE)

rm(static_var_top)
