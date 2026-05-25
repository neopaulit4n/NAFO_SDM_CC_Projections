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
  )

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

# SST Climate Projection Table — generalised for any number of variables
# Input dataframe columns:
#   period     : "Reference", "P1", "P2", "P3", "P4"
#   ssp        : "Reference", "SSP 1-2.6", "SSP 2-4.5", "SSP 3-7.0", "SSP 5-8.5"
#   variable   : e.g. "SST range", "SST min", …  (becomes a row-group header)
#   mean       : numeric
#   sd         : numeric
#
# Optional extra columns (kept as-is, appended after the period columns):
#   adf        : numeric  — ADF test statistic
#   conclusion : character — "Stationary" / "Non-stationary"
#
# Install if needed:
#   install.packages(c("gt", "dplyr", "tidyr"))

# ── 0.  Customise these mappings ──────────────────────────────────────────────

# Human-readable column labels for each period value
PERIOD_LABELS <- c(
  "Reference" = "1993\u20132014", # en-dash via unicode
  "P1" = "P1: 2020\u20132039",
  "P2" = "P2: 2040\u20132059",
  "P3" = "P3: 2060\u20132079",
  "P4" = "P4: 2080\u20132099"
)

# Desired column order for the period columns in the final table
PERIOD_ORDER <- c("Reference", "P1", "P2", "P3", "P4")

# Number of decimal places for mean ± sd cells
N_DECIMALS <- 2

# ── 1.  Example input data (replace with your real data frame) ────────────────

df_input <- tibble::tribble(
  ~period,      ~ssp,          ~variable,    ~mean,  ~sd,   ~adf,   ~conclusion,
  "Reference",  "Reference",   "SST range",  12.74,  2.08,   NA,    NA,
  "P1",         "SSP 1-2.6",   "SST range",  12.85,  1.91,  -1.12,  "Non-stationary",
  "P2",         "SSP 1-2.6",   "SST range",  12.99,  1.83,  -1.12,  "Non-stationary",
  "P3",         "SSP 1-2.6",   "SST range",  13.13,  1.67,  -1.12,  "Non-stationary",
  "P4",         "SSP 1-2.6",   "SST range",  13.09,  1.70,  -1.12,  "Non-stationary",
  "P1",         "SSP 2-4.5",   "SST range",  13.12,  1.97,  -2.12,  "Non-stationary",
  "P2",         "SSP 2-4.5",   "SST range",  13.36,  1.93,  -2.12,  "Non-stationary",
  "P3",         "SSP 2-4.5",   "SST range",  13.29,  1.79,  -2.12,  "Non-stationary",
  "P4",         "SSP 2-4.5",   "SST range",  13.69,  1.67,  -2.12,  "Non-stationary",
  "P1",         "SSP 3-7.0",   "SST range",  12.59,  1.87,  -9.89,  "Stationary",
  "P2",         "SSP 3-7.0",   "SST range",  13.35,  2.00,  -9.89,  "Stationary",
  "P3",         "SSP 3-7.0",   "SST range",  13.85,  1.92,  -9.89,  "Stationary",
  "P4",         "SSP 3-7.0",   "SST range",  14.53,  1.84,  -9.89,  "Stationary",
  "P1",         "SSP 5-8.5",   "SST range",  13.07,  2.06,   1.24,  "Non-stationary",
  "P2",         "SSP 5-8.5",   "SST range",  13.69,  1.80,   1.24,  "Non-stationary",
  "P3",         "SSP 5-8.5",   "SST range",  14.47,  1.72,   1.24,  "Non-stationary",
  "P4",         "SSP 5-8.5",   "SST range",  16.24,  1.69,   1.24,  "Non-stationary",

  "Reference",  "Reference",   "SST min",    1.57,   1.40,   NA,    NA,
  "P1",         "SSP 1-2.6",   "SST min",    2.47,   1.47,  -6.41,  "Stationary",
  "P2",         "SSP 1-2.6",   "SST min",    2.72,   1.42,  -6.41,  "Stationary",
  "P3",         "SSP 1-2.6",   "SST min",    2.73,   1.40,  -6.41,  "Stationary",
  "P4",         "SSP 1-2.6",   "SST min",    2.83,   1.38,  -6.41,  "Stationary",
  "P1",         "SSP 2-4.5",   "SST min",    2.33,   1.43,  -1.77,  "Non-stationary",
  "P2",         "SSP 2-4.5",   "SST min",    2.71,   1.43,  -1.77,  "Non-stationary",
  "P3",         "SSP 2-4.5",   "SST min",    3.18,   1.39,  -1.77,  "Non-stationary",
  "P4",         "SSP 2-4.5",   "SST min",    3.43,   1.40,  -1.77,  "Non-stationary",
  "P1",         "SSP 3-7.0",   "SST min",    2.45,   1.40,  -7.67,  "Stationary",
  "P2",         "SSP 3-7.0",   "SST min",    2.85,   1.41,  -7.67,  "Stationary",
  "P3",         "SSP 3-7.0",   "SST min",    3.39,   1.41,  -7.67,  "Stationary",
  "P4",         "SSP 3-7.0",   "SST min",    3.97,   1.43,  -7.67,  "Stationary",
  "P1",         "SSP 5-8.5",   "SST min",    2.40,   1.42,  -3.52,  "Non-stationary",
  "P2",         "SSP 5-8.5",   "SST min",    2.92,   1.38,  -3.52,  "Non-stationary",
  "P3",         "SSP 5-8.5",   "SST min",    3.77,   1.34,  -3.52,  "Non-stationary",
  "P4",         "SSP 5-8.5",   "SST min",    4.42,   1.37,  -3.52,  "Non-stationary"
)

# ── 2.  Prepare: tidy long → wide ────────────────────────────────────────────

df <- df4 |>
  mutate(
    cell = sprintf(paste0("%.", 2, "f \u00b1 %.", 2, "f"), mean, sd),
    cell = ifelse(bold, paste0("<b>",cell,"</b>"), cell),
    ADF = ifelse(period == "Reference" & ssp == "Reference", NA, ADF),
    Conclusion = ifelse(period == "Reference" & ssp == "Reference", NA, Conclusion),
  ) |>
  select(-mean, -sd, -bold) |>
  distinct() |>
  pivot_wider(names_from = period, values_from = cell) |>
  relocate(ADF:Conclusion, .after = last_col())

prepare_gt_data <- function(df, period_order = PERIOD_ORDER, n_dec = N_DECIMALS) {

  # Detect which optional extra columns exist (anything beyond the core 5)
  core_cols  <- c("period", "ssp", "variable", "mean", "sd")
  extra_cols <- setdiff(names(df), core_cols)   # e.g. c("adf", "conclusion")

  # Build the mean ± sd cell string
  df <- df |>
    mutate(
      cell = sprintf(paste0("%.", n_dec, "f \u00b1 %.", n_dec, "f"), mean, sd)
    )

  # --- Pivot period columns wide -------------------------------------------
  wide_periods <- df |>
    select(variable, ssp, period, cell) |>
    pivot_wider(names_from = period, values_from = cell) |>
    # Ensure only the periods we expect are present (and in order)
    select(variable, ssp, any_of(period_order))

  # --- Collapse extra cols: one row per (variable × ssp) -------------------
  # For numeric extras (e.g. adf) take the first non-NA value;
  # for character extras (e.g. conclusion) take the first non-NA value.
  if (length(extra_cols) > 0) {
    wide_extras <- df |>
      group_by(variable, ssp) |>
      summarise(
        across(all_of(extra_cols),
               ~ if (is.numeric(.x)) first(na.omit(.x)) else first(na.omit(.x))),
        .groups = "drop"
      )
    wide <- left_join(wide_periods, wide_extras, by = c("variable", "ssp"))
  } else {
    wide <- wide_periods
  }

  # --- Reorder SSP rows: Reference first, then SSP rows in their order ------
  ssp_levels <- c("Reference", setdiff(unique(df$ssp), "Reference"))
  wide <- wide |>
    mutate(ssp = factor(ssp, levels = ssp_levels)) |>
    arrange(variable, ssp) |>
    mutate(ssp = as.character(ssp))

  # --- Rename ssp → SSP for display ----------------------------------------
  wide <- rename(wide, SSP = ssp)

  wide
}

gt_data <- prepare_gt_data(df_input)

# ── 3.  Detect which period and extra columns are actually present ─────────────

present_periods <- intersect(PERIOD_ORDER, names(gt_data))
extra_cols      <- setdiff(names(gt_data), c("variable", "SSP", present_periods))
has_adf        <- "adf"        %in% extra_cols
has_conclusion <- "conclusion" %in% extra_cols

# ── 4.  Build the gt table ────────────────────────────────────────────────────

build_gt <- function(data,
                     present_periods,
                     period_labels  = PERIOD_LABELS,
                     extra_cols     = character(0)) {

  # Dynamic column-label list for period columns
  period_label_list <- setNames(
    as.list(unname(period_labels[present_periods])),
    present_periods
  )

  # Dynamic column-label list for extra cols (capitalise first letter)
  extra_label_list <- setNames(
    as.list(tools::toTitleCase(extra_cols)),
    extra_cols
  )

  all_label_list <- c(list(SSP = "SSP"), period_label_list, extra_label_list)

  tbl <- data |>
    gt(groupname_col = "variable") |>

    # ── Column labels ────────────────────────────────────────────────────────
    cols_label(.list = all_label_list) |>

    # ── Spanner over period columns ──────────────────────────────────────────
    tab_spanner(
      label   = "Time Period",
      columns = all_of(present_periods)
    ) |>

    # ── Replace NAs ─────────────────────────────────────────────────────────
    sub_missing(missing_text = "") |>

    # ── Alignment ───────────────────────────────────────────────────────────
    cols_align(align = "center", columns = all_of(present_periods)) |>
    cols_align(align = "left",   columns = SSP)

  # Optional: centre ADF, left-align Conclusion
  if (has_adf)        tbl <- tbl |> cols_align(align = "center", columns = adf)
  if (has_conclusion) tbl <- tbl |> cols_align(align = "left",   columns = conclusion)

  # ── Styling ─────────────────────────────────────────────────────────────────
  tbl <- tbl |>

    tab_style(
      style = list(cell_fill(color = "#f0f0f0"),
                   cell_text(weight = "bold", style = "italic")),
      locations = cells_row_groups()
    ) |>

    tab_style(
      style = cell_text(weight = "bold"),
      locations = cells_column_labels()
    ) |>

    tab_style(
      style = cell_text(weight = "bold"),
      locations = cells_column_spanners()
    ) |>

    tab_style(
      style = cell_fill(color = "#e8f4fd"),
      locations = cells_body(rows = SSP == "Reference")
    )

  # Conditional conclusion colouring (only if column exists)
  if (has_conclusion) {
    tbl <- tbl |>
      tab_style(
        style = cell_text(color = "#2e7d32", weight = "bold"),
        locations = cells_body(columns = conclusion,
                               rows    = conclusion == "Stationary")
      ) |>
      tab_style(
        style = cell_text(color = "#b45309"),
        locations = cells_body(columns = conclusion,
                               rows    = conclusion == "Non-stationary")
      )
  }

  # ── Table options ────────────────────────────────────────────────────────────
  tbl <- tbl |>
    tab_options(
      table.font.size                    = px(13),
      column_labels.font.size            = px(13),
      row_group.font.size                = px(13),
      heading.align                      = "left",
      table.border.top.style             = "solid",
      table.border.top.width             = px(2),
      table.border.top.color             = "#333333",
      table.border.bottom.style          = "solid",
      table.border.bottom.width          = px(2),
      table.border.bottom.color          = "#333333",
      column_labels.border.top.style     = "solid",
      column_labels.border.top.width     = px(1),
      column_labels.border.top.color     = "#999999",
      column_labels.border.bottom.style  = "solid",
      column_labels.border.bottom.width  = px(2),
      column_labels.border.bottom.color  = "#333333",
      row_group.border.top.style         = "solid",
      row_group.border.top.width         = px(1),
      row_group.border.top.color         = "#cccccc",
      row_group.border.bottom.style      = "hidden",
      stub.border.style                  = "hidden",
      data_row.padding                   = px(5),
      table.width                        = pct(100)
    ) |>

    tab_caption(
      caption = md("**Table X.** Projected sea surface temperature statistics by SSP scenario and time period. ADF = Augmented Dickey-Fuller test statistic.")
    )

  tbl
}

tbl <- build_gt(gt_data, present_periods, extra_cols = extra_cols)

# ── 5.  Save ─────────────────────────────────────────────────────────────────

gtsave(tbl, filename = "sst_table.html")

# PNG (requires webshot2):
# gtsave(tbl, filename = "sst_table.png", zoom = 2)

message("\u2713  Table saved to sst_table.html")