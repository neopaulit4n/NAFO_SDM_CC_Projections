# Create pretty table replicating table 5 from black corals
df1 <- bind_cols(z,zz) |>
  pivot_longer(
    cols = matches("^P\\d"),
    names_to = c("period", "ssp", ".value"),
    names_pattern = "^(P\\d)\\.(\\d-\\d+\\.\\d)\\.(.+)$"
  ) |>
  mutate(
    period = ifelse(period == "P0", "Reference", period),
    ssp = ifelse(ssp == "0-0.0", "Reference", ssp)
  ) |>
  pivot_longer(cols = -c("period","ssp"), names_to = "variable", values_to = "value")

df2 <- lapply(selected_cmip_vars, function(var) {
  lapply(ssp_all, function(sspoi) {
    df <- filter(df1, ssp == sspoi, variable == var)
    df_aov <- aov(value ~ period, data = df)
    df_tukey <- TukeyHSD(df_aov, conf.level = 0.95)
    df <- data.frame(
      period = period_all,
      bold = c(
        ifelse(df_tukey$period[1,4] > 0.05, TRUE, FALSE),  # P1
        ifelse(df_tukey$period[1,4] > 0.05 | df_tukey$period[4,4] > 0.05, TRUE, FALSE),  # P2
        ifelse(df_tukey$period[4,4] > 0.05 | df_tukey$period[6,4] > 0.05, TRUE, FALSE),  # P3
        ifelse(df_tukey$period[6,4] > 0.05, TRUE, FALSE)  # P4
      ),
      ssp = sspoi,
      variable = var
    )
  }) |>
    bind_rows()
}) |>
  bind_rows()

df3 <- left_join(df1, df2, by = c("period", "ssp", "variable")) |>
  summarise(
    mean = mean(value, na.rm = TRUE), 
    sd = sd(value, na.rm = TRUE), 
    .by = c(period, ssp, variable, bold)
  ) |>
  mutate(bold = ifelse(is.na(bold), FALSE, bold))

df4 <- lapply(selected_cmip_vars, function(var) {
  lapply(ssp_all, function(sspoi) {
    df <- filter(df3, ssp %in% c("Reference", sspoi), variable == var)
    adf_res <- tseries::adf.test(df$mean, k = 0)
    df$ADF = adf_res$statistic
    df$Conclusion = ifelse(adf_res$p.value < 0.05, "Non-stationary", "Stationary")
    return(df)
  }) |>
    bind_rows()
}) |>
  bind_rows() #|>
  # mutate(ssp = ifelse(ssp == "Reference", NA, ssp)) |>
  # fill(ssp, .direction = "up")



df <- df4 |>
  mutate(
    cell = sprintf(paste0("%.", 2, "f \u00b1 %.", 2, "f"), mean, sd),
    cell = ifelse(bold, paste0("<b>",cell,"</b>"), cell),
    ADF = ifelse(period == "Reference" & ssp == "Reference", NA, round(ADF, 2)),
    Conclusion = ifelse(period == "Reference" & ssp == "Reference", NA, Conclusion),
    period = str_replace_all(period, c(
      "Reference" = "1993-2014",
      "P1" = "P1: 2020-2039",
      "P2" = "P2: 2040-2059",
      "P3" = "P3: 2060-2079",
      "P4" = "P4: 2080-2099"
    ))
  ) |>
  select(-mean, -sd, -bold) |>
  distinct() |>
  pivot_wider(names_from = period, values_from = cell) |>
  relocate(ADF:Conclusion, .after = last_col()) |>
  mutate(
    var_clean = str_replace_all(variable, c(
      # Variables
      "^BS" = "Bottom Salinity",
      "^SSS" = "Sea Surface Salinity",
      "^BT" = "Bottom Temperature",
      "^SST" = "Sea Surface Temperature",
      "^BCS" = "Bottom Current Speed",
      "^BStr" = "Bottom Stress",
      "^MLD_W" = "Winter Mixed Layer Depth",
      "^MLD_Sp" = "Spring Mixed Layer Depth",
      "^MLD_Su" = "Summer Mixed Layer Depth",
      "^MLD_F" = "Fall Mixed Layer Depth",
      # "^MLD" = "Annual Mixed Layer Depth"  # abbreviation will need to be changed to "MLD_An"
      
      # Statistics
      "_min" = " Minimum",
      "_max" = " Maximum",
      "_mean" = " Mean",
      "_range" = " Range"
    )),
    var_abbr = str_replace_all(variable, "_", " "),
    var_abbr = paste0("(",var_abbr,")"),
    var_suffix = case_when(
      str_detect(var_clean, "Salinity") ~ "",
      str_detect(var_clean, "Temperature") ~ "(°C)",
      str_detect(var_clean, "Speed") ~ "(m/s)",
      str_detect(var_clean, "Stress") ~ "(PA)",
      str_detect(var_clean, "Depth") ~ "(m)"
    ),
    var_final = paste(var_clean, var_abbr, var_suffix)
  ) |>
    select(-c(variable, var_clean, var_abbr, var_suffix)) |>
    relocate(var_final, .after = ssp) |>
    rename(SSP = ssp)

tbl <- df |>
  gt::gt(groupname_col = "var_final") |>
  gt::tab_spanner(
    label   = "Time Period",
    columns = c("1993-2014","P1: 2020-2039","P2: 2040-2059","P3: 2060-2079","P4: 2080-2099")
  ) |>
  gt::sub_missing(missing_text = "") |>
  gt::cols_align(align = "center", columns = c("1993-2014","P1: 2020-2039","P2: 2040-2059","P3: 2060-2079","P4: 2080-2099")) |>
  gt::cols_align(align = "left",   columns = SSP) |>
  gt::fmt_markdown() |>  # applying bold
  gt::tab_options(
    table_body.hlines.style = "none",
    table_body.vlines.style = "none"
  ) |>
  flextable::as_flextable() |>
  flextable::border_remove() |>
  flextable::save_as_docx(path = paste0(output_folder,"/",vmeoi,"_ADFTukeyCMIPVarTable2.docx"))

gt::gtsave(tbl, paste0(output_folder,"/",vmeoi,"_ADFTukeyCMIPVarTable.docx"))

doc <- officer::read_docx(paste0(output_folder,"/",vmeoi,"_ADFTukeyCMIPVarTable.docx"))
doc_xml <- officer::docx_body_xml(doc)
doc_bord <- xml2::xml_find_all(doc_xml, ".//w:tcBorders")

for (node in doc_bord) {
  tops    <- xml2::xml_find_all(node, "w:top")
  bottoms <- xml2::xml_find_all(node, "w:bottom")
  starts  <- xml2::xml_find_all(node, "w:start")
  ends    <- xml2::xml_find_all(node, "w:end")
  
  for (el in c(tops, bottoms, starts, ends)) {
    sz <- xml2::xml_attr(el, "sz")
    if (is.na(sz)) {
      xml2::xml_remove(el)
    }
  }
}

print(doc, target = paste0(output_folder,"/",vmeoi,"_ADFTukeyCMIPVarTable.docx"))
