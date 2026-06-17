# Create pretty tables of extrapolation area percentages

# Read in extrapolation dataframes ----
ea_files <- list.files(
  path = paste0(output_folder,"/Extrapolations"),
  pattern = "csv$",
  full.names = TRUE)

ea <- lapply(ea_files, read_csv, show_col_types = FALSE) |>
  set_names(gsub("^.+?_(.+?)\\.csv$", "\\1", ea_files)) |>
  bind_rows(.id = "PeriodSSP")

rm(ea_files)

# Table 7: Percentages of areas of extrapolation ----
ea_overall <- filter(ea, covariate == "Overall") |>
  select(-freq, -covariate) |>
  mutate(
    Period = ifelse(PeriodSSP == "baseline", "Reference", gsub("(P\\d)\\.\\d-\\d\\.\\d", "\\1", PeriodSSP)),
    SSP = ifelse(PeriodSSP == "baseline", "", gsub("P\\d\\.(\\d-\\d\\.\\d)", "\\1", PeriodSSP))
  ) |>
  select(-PeriodSSP) |>
  mutate(
    InputData = replace_values(InputData, from = c("PA", "PresenceOnly"), to = c("Presence and Absence Data", "Presence Data")),
    Type = paste("Percent", Type, "Extrapolation")
  ) |>
  pivot_wider(id_cols = c("Period","SSP"), names_from = c("Type", "InputData"), values_from = "perc") |>
  mutate(
    across(where(is.double), ~ ifelse(is.na(.x), 0, .x)),
    Period = ifelse(duplicated(Period), "", Period))



ea_tbl7 <- gt::gt(ea_overall) |>
  gt::tab_spanner(
    label = "Presence and Absence Data",
    columns = ends_with("Absence Data"),
    id = "spanner1") |>
  gt::tab_spanner(
    label = "Presence Data",
    columns = ends_with("Presence Data"),
    id = "spanner2") |>
  gt::cols_label(
    ends_with("Absence Data") ~ "",
    ends_with("Presence Data") ~ ""
  )

gt::gtsave(ea_tbl7, paste0(output_folder,"/Extrapolations/",vmeoi,"_OverallPercentagesTable.docx"))
