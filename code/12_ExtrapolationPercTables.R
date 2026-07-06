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
ea_tbl7 <- filter(ea, covariate == "Overall") |>
  select(-freq, -covariate) |>
  mutate(
    Period = ifelse(PeriodSSP == "baseline", "Reference", gsub("(P\\d)\\.\\d-\\d\\.\\d", "\\1", PeriodSSP)),
    SSP = ifelse(PeriodSSP == "baseline", "", gsub("P\\d\\.(\\d-\\d\\.\\d)", "\\1", PeriodSSP))) |>
  select(-PeriodSSP) |>
  mutate(
    InputData = replace_values(InputData, from = c("PA", "PresenceOnly"), to = c("Presence and Absence Data", "Presence Data")),
    Type = paste("Percent", Type, "Extrapolation"),
    perc_modif = ifelse(perc < 0.01, "less_than", "round"),
    perc = ifelse(perc_modif == "less_than", "< 0.01", as.character(round(perc, 2)))) |>
  pivot_wider(id_cols = c("Period","SSP"), names_from = c("Type", "InputData"), values_from = "perc") |>
  mutate(
    across(starts_with("Percent"), ~ ifelse(is.na(.x), "0", .x)),
    Period = ifelse(duplicated(Period), "", Period)) |>
  gt::gt() |>
  gt::tab_spanner(
    label = "Presence and Absence Data",
    columns = ends_with("Absence Data"),
    id = "spanner1") |>
  gt::tab_spanner(
    label = "Presence Data",
    columns = ends_with("Presence Data"),
    id = "spanner2") |>
  gt::cols_label(
    "Percent Univariate Extrapolation_Presence and Absence Data" ~ "Percent Univariate Extrapolation",
    "Percent Combinatorial Extrapolation_Presence and Absence Data" ~ "Percent Combinatorial Extrapolation",
    "Percent Analogue Extrapolation_Presence and Absence Data" ~ "Percent Analogue Extrapolation",
    "Percent Univariate Extrapolation_Presence Data" ~ "Percent Univariate Extrapolation",
    "Percent Combinatorial Extrapolation_Presence Data" ~ "Percent Combinatorial Extrapolation",
    "Percent Analogue Extrapolation_Presence Data" ~ "Percent Analogue Extrapolation") |>
  gt::tab_style(
    style = gt::cell_text(weight = "bold"),
    locations = list(gt::cells_column_spanners(), gt::cells_column_labels())
  ) |>
  gt::cols_align(align = "left", columns = starts_with("Percent"))

ea_tbl7

if (file.exists(paste0(output_folder,"/Extrapolations/",vmeoi,"_OverallPercentagesTable.docx"))) file.remove(paste0(output_folder,"/Extrapolations/",vmeoi,"_OverallPercentagesTable.docx"))
gt::gtsave(ea_tbl7, paste0(output_folder,"/Extrapolations/",vmeoi,"_OverallPercentagesTable.docx"))

doc <- officer::read_docx(paste0(output_folder, "/Extrapolations/", vmeoi, "_OverallPercentagesTable.docx"))
doc_xml <- officer::docx_body_xml(doc)

# Get all rows in the table
all_rows <- xml2::xml_find_all(doc_xml, ".//w:tr")

# Identify which row indices contain period labels
period_labels <- c("Reference", "P1", "P2", "P3", "P4")

period_row_indices <- c()
for (i in seq_along(all_rows)) {
  row_text <- xml2::xml_text(all_rows[[i]])
  if (any(trimws(row_text) == period_labels | grepl(paste(period_labels, collapse = "|"), row_text))) {
    period_row_indices <- c(period_row_indices, i)
  }
}

# Now process all cell borders
doc_bord <- xml2::xml_find_all(doc_xml, ".//w:tcBorders")

for (node in doc_bord) {
  # Find which row this border node belongs to
  parent_row <- xml2::xml_find_first(node, "./ancestor::w:tr")
  row_index  <- which(sapply(all_rows, identical, parent_row))
  in_period_row <- length(row_index) > 0 && row_index %in% period_row_indices

  tops    <- xml2::xml_find_all(node, "w:top")
  bottoms <- xml2::xml_find_all(node, "w:bottom")
  starts  <- xml2::xml_find_all(node, "w:start")
  ends    <- xml2::xml_find_all(node, "w:end")

  for (el in bottoms) {
    sz <- xml2::xml_attr(el, "sz")
    if (is.na(sz)) xml2::xml_remove(el)
  }
  for (el in c(starts, ends)) {
    sz <- xml2::xml_attr(el, "sz")
    if (is.na(sz)) xml2::xml_remove(el)
  }
  for (el in tops) {
    sz <- xml2::xml_attr(el, "sz")
    if (is.na(sz) && !in_period_row) {
      # Remove top border only if NOT a period row
      xml2::xml_remove(el)
    }
  }
}

print(doc, target = paste0(output_folder, "/Extrapolations/", vmeoi, "_OverallPercentagesTable.docx"))

# Supplementary table S3 ----
ea_tbls3 <- filter(ea, covariate != "Overall", Type == "Univariate") |>
  select(-freq, -Type) |>
  complete(PeriodSSP, InputData, covariate) |>
  mutate(
    Period = ifelse(PeriodSSP == "baseline", "Reference", gsub("(P\\d)\\.\\d-\\d\\.\\d", "\\1", PeriodSSP)),
    Period = factor(Period, levels = period_labels),
    SSP = ifelse(PeriodSSP == "baseline", "", gsub("P\\d\\.(\\d-\\d\\.\\d)", "\\1", PeriodSSP))#,
    # perc = ifelse(is.na(perc), 0, perc)
  ) |>
  select(-PeriodSSP) |>
  group_by(Period, SSP, InputData) |>
  arrange(desc(perc), .by_group = TRUE) |>
  mutate(rank = row_number()) |>
  ungroup() |>
  mutate(
    InputData = replace_values(InputData, from = c("PA", "PresenceOnly"), to = c("Presence and Absence Data", "Presence Data")),
    perc_modif = ifelse(perc < 0.01, "less_than", "round"),
    perc = ifelse(perc_modif == "less_than", "< 0.01", as.character(round(perc, 2)))
  ) |>
  pivot_wider(id_cols = c(Period, SSP, rank), names_from = InputData, values_from = c(covariate, perc), names_glue = "{InputData}_{.value}") |>
  select(-rank) |>
  mutate(across(starts_with("Percent"), ~ ifelse(is.na(.x), "0", .x))
  ) |>
  mutate(SSP = ifelse(duplicated(SSP), "", SSP), .by = Period) |>
  mutate(Period = ifelse(duplicated(Period), "", as.character(Period))) |>
  relocate(`Presence Data_covariate`, .before = `Presence Data_perc`) |>
  gt::gt() |>
  gt::tab_spanner(
    label = "Presence and Absence Data",
    columns = starts_with("Presence and Absence"),
    id = "spanner1") |>
  gt::tab_spanner(
    label = "Presence Data",
    columns = starts_with("Presence Data"),
    id = "spanner2") |>
  gt::cols_label(
    ends_with("covariate") ~ "Variable",
    ends_with("perc") ~ "Percent Area"
  ) |>
  gt::tab_style(
    style = gt::cell_text(weight = "bold"),
    locations = list(gt::cells_column_spanners(), gt::cells_column_labels())
  ) |>
  gt::cols_align(align = "left", columns = starts_with("Percent"))

ea_tbls3

if (file.exists(paste0(output_folder,"/Extrapolations/",vmeoi,"_VariablePercentagesTable.docx"))) file.remove(paste0(output_folder,"/Extrapolations/",vmeoi,"_VariablePercentagesTable.docx"))
gt::gtsave(ea_tbls3, paste0(output_folder,"/Extrapolations/",vmeoi,"_VariablePercentagesTable.docx"))

doc <- officer::read_docx(paste0(output_folder, "/Extrapolations/", vmeoi, "_VariablePercentagesTable.docx"))
doc_xml <- officer::docx_body_xml(doc)

# Get all rows in the table
all_rows <- xml2::xml_find_all(doc_xml, ".//w:tr")

# Identify which row indices contain period labels
period_labels <- c("Reference", "P1", "P2", "P3", "P4", "1-2.6", "2-4.5", "3-7.0", "5-8.5")

period_row_indices <- c()
for (i in seq_along(all_rows)) {
  row_text <- xml2::xml_text(all_rows[[i]])
  if (any(trimws(row_text) == period_labels | grepl(paste(period_labels, collapse = "|"), row_text))) {
    period_row_indices <- c(period_row_indices, i)
  }
}

# Now process all cell borders
doc_bord <- xml2::xml_find_all(doc_xml, ".//w:tcBorders")

for (node in doc_bord) {
  # Find which row this border node belongs to
  parent_row <- xml2::xml_find_first(node, "./ancestor::w:tr")
  row_index  <- which(sapply(all_rows, identical, parent_row))
  in_period_row <- length(row_index) > 0 && row_index %in% period_row_indices

  tops    <- xml2::xml_find_all(node, "w:top")
  bottoms <- xml2::xml_find_all(node, "w:bottom")
  starts  <- xml2::xml_find_all(node, "w:start")
  ends    <- xml2::xml_find_all(node, "w:end")

  for (el in bottoms) {
    sz <- xml2::xml_attr(el, "sz")
    if (is.na(sz)) xml2::xml_remove(el)
  }
  for (el in c(starts, ends)) {
    sz <- xml2::xml_attr(el, "sz")
    if (is.na(sz)) xml2::xml_remove(el)
  }
  for (el in tops) {
    sz <- xml2::xml_attr(el, "sz")
    if (is.na(sz) && !in_period_row) {
      # Remove top border only if NOT a period row
      xml2::xml_remove(el)
    }
  }
}
print(doc, target = paste0(output_folder, "/Extrapolations/", vmeoi, "_VariablePercentagesTable.docx"))