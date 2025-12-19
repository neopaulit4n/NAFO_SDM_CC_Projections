
# Add VME_P_A field to response dataframes that are missing it

# Load raw VME dataframes
resp_list <- list.files(path = "data/raw/SDM2024/FINAL Response Variables",
                      pattern = "*.csv", recursive = TRUE, full.names = TRUE) %>%
  set_names(., nm = basename(.) %>% tools::file_path_sans_ext()) %>%
  lapply(read_csv, show_col_types = FALSE) %>%
  lapply(function(df) {
    df <- mutate(df, Set = as.character(Set))
    return(df)
  })

# Check out the dataframes to see what they have in common/what needs to be changed
sapply(resp_list, print(head))

# Coalesce dataframes
resp_df <- bind_rows(resp_list, .id = "VME_Group") %>%
  mutate(VME_Group_kg = coalesce(LargeGorgoniansFG,SeaPensFG,SmallGorgoniansFG,LargeSpongesFG_v2,SpongesFG),
         VME_Group = case_match(VME_Group,
                                "boltenia_sdm_response" ~ "boltenia",
                                "bryozoan_sdm_response" ~ "bryozoan",
                                "lg_sdm_response" ~ "large_gorgonians",
                                "sg_sdm_response" ~ "small_gorgonians",
                                "largespongesFG" ~ "large_sponges",
                                .default = VME_Group
                                )
         ) %>%
  select(-c(ID,VME_Group_2,VME_threshold,
            LargeGorgoniansFG,SeaPensFG,SmallGorgoniansFG,LargeSpongesFG_v2,SpongesFG)) %>%
  rename(Other_seapens = Other_taxa) %>%
  pivot_longer(cols = c(Anthoptilum,Balticina,Funiculina,Pennatula,Other_seapens,Acanella,
                        Anthothelidae,Radicipes,Astrophorina,Tetillidae,Polymastiidae),
               names_to = "Taxa2",
               values_to = "Taxa_Biomass_Kg2") %>%
  mutate(Taxa = coalesce(Taxa, Taxa2),
         Taxa_Biomass_Kg = coalesce(Taxa_Biomass_Kg, Taxa_Biomass_Kg2)) %>%
  select(-c(Taxa2, Taxa_Biomass_Kg2)) %>%
  mutate(Taxa = ifelse(is.na(Taxa_Biomass_Kg), NA, Taxa)) %>%
  distinct %>%
  group_by(VME_Group, Mission, Set, Start_Lat_DD, Start_Long_DD) %>%
  mutate(VME_Group_kg2 = sum(Taxa_Biomass_Kg, na.rm = TRUE),
         VME_Group_kg = coalesce(VME_Group_kg, VME_Group_kg2),
         VME_P_A2 = ifelse(VME_Group_kg > 0, 1, 0),
         VME_P_A = coalesce(VME_P_A, VME_P_A2)) %>%
  ungroup %>%
  select(-c(Taxa, Taxa_Biomass_Kg, VME_Group_kg2, VME_P_A2)) %>%
  distinct
write_csv(resp_df, "data/cleaned/VME_group_PA_df.csv")

# If wanted to keep taxa, would need to resolve boltenia situation where some rows don't have any taxa/kg data
  # but still have VME_P_A = 1. Some taxa for boltenia also marked as TUNICATE, SESSILE - does that still count as boltenia?
  # Use group_by(VME_Group) %>% complete(nesting(all_the_set_information), Taxa, fill = list(Taxa_Biomass_Kg = 0)) %>%
  # ungroup() to fill in missing combinations with 0 biomass




