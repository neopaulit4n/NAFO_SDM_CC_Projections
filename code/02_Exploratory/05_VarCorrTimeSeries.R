
# Exploring how variable correlations change across time periods (stationarity of variable correlations)

# Create directories ----
if (!dir.exists(paste0(output_folder,"/VarCorrelations"))) dir.create(paste0(output_folder,"/VarCorrelations"))
if (!dir.exists(paste0(output_folder,"/VarCorrelations/01_InputDataframes"))) dir.create(paste0(output_folder,"/VarCorrelations/01_InputDataframes"))
if (!dir.exists(paste0(output_folder,"/VarCorrelations/02_CorMat"))) dir.create(paste0(output_folder,"/VarCorrelations/02_CorMat"))
if (!dir.exists(paste0(output_folder,"/VarCorrelations/03_CorPlots"))) dir.create(paste0(output_folder,"/VarCorrelations/03_CorPlots"))
if (!dir.exists(paste0(output_folder,"/VarCorrelations/04_CorDiffMat"))) dir.create(paste0(output_folder,"/VarCorrelations/04_CorDiffMat"))
if (!dir.exists(paste0(output_folder,"/VarCorrelations/05_CorDiffPlots"))) dir.create(paste0(output_folder,"/VarCorrelations/05_CorDiffPlots"))

# Plot correlations between variable pairs over time
baseline_layers <- c(cmip_layers, terrain = bathy_layers[[vme_terrain_vars]]) |>
  terra::rast()
names(baseline_layers) <- c(names(cmip_layers), vme_terrain_vars)

projection_layers <- unlist(cmip_layers_proj, recursive = FALSE) |>
  lapply(X = _, 
    function(layer) {
      rast_result <- c(layer, bathy_layers[[vme_terrain_vars]]) |>
        terra::rast()
      names(rast_result) <- gsub("GEBCO2024_FS005_StudyArea_","", names(rast_result))
      return(rast_result)
    })

all_layers <- list()
all_layers[[1]] <- baseline_layers
names(all_layers) <- "baseline"
all_layers <- c(all_layers, projection_layers)
all_layers_df <- lapply(all_layers, function(x) terra::as.data.frame(x) |> drop_na())

# Output base dataframes used to calculate correlations between variables
lapply(1:length(all_layers_df), function(df) {
  write_csv(all_layers_df[[df]], paste0(output_folder,"/VarCorrelations/01_InputDataframes/FullArea_",names(all_layers_df)[df],"_df.csv"))
})

# Calculate correlations ----
all_layers_cor <- lapply(1:length(all_layers_df), function(layer) {
  cor_mat <- cor(all_layers_df[[layer]], method = "spearman", use = "complete.obs") |>
    as.data.frame()
  write.csv(cor_mat, paste0(output_folder,"/VarCorrelations/02_CorMat/FullArea_",names(all_layers_df)[layer],"_cor.csv"))
  cor_long <- as.data.frame(cor_mat) %>%
    rownames_to_column(var = "var1") %>%
    pivot_longer(-var1, names_to = "var2", values_to = "cor")
  return(cor_long)
}) %>%
  set_names(names(all_layers_df)) %>%
  bind_rows(.id = "PeriodSSP") %>%
  mutate(
    period = ifelse(is.na(str_extract(PeriodSSP, "P\\d")), "baseline", str_extract(PeriodSSP, "P\\d")),
    ssp = ifelse(is.na(str_extract(PeriodSSP, "\\d-\\d\\.\\d")), "", str_extract(PeriodSSP, "\\d-\\d\\.\\d"))
  )

# Output correlation plots ----
lapply(names(all_layers_df), function(x) {
  df <- filter(all_layers_cor, PeriodSSP == x)
  x_name <- ifelse(x == "baseline", "Reference", paste(df$period[1], "SSP", df$ssp[1]))
  p <- ggplot(data = df, aes(x = var1, y = var2, fill = cor)) +
    geom_tile(colour = "black") +
    geom_text(aes(label = ifelse(cor > 0.7 | cor < -0.7, "*", ""))) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", limit = c(-1,1), midpoint = 0) +
    scale_x_discrete(expand = c(0,0)) +
    scale_y_discrete(expand = c(0,0)) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text = element_text(size = 8),
          axis.line = element_blank(),
          axis.title = element_blank()) +
    labs(fill = "Correlation",
        title = x_name)
  ggsave(
    filename = paste0(output_folder,"/VarCorrelations/03_CorPlots/FullArea_",x,".jpg"),
    plot = p, width = 8, height = 6
  )
})

# Do the same but only for cells with overlapping VME P/A data ----
all_layers_overlap <- lapply(all_layers, function(layer) {
  terra::extract(
    layer, 
    select(
      filter(resp_df, VME_Group == vmeoi), #%>%
        #drop_na(),
      Start_Long_DD, Start_Lat_DD)
    ) %>%
    select(-ID)
}) 

# Output base dataframes used to calculate correlations between variables
lapply(1:length(all_layers_overlap), function(df) {
  write_csv(all_layers_overlap[[df]], paste0(output_folder,"/VarCorrelations/01_InputDataframes/OverlapVME_",names(all_layers_overlap)[df],"_df.csv"))
})

# Calculate correlations ----
all_layers_overlap_cor <- lapply(1:length(all_layers_overlap), function(layer) {
  cor_mat <- cor(all_layers_overlap[[layer]], method = "spearman", use = "complete.obs") |>
    as.data.frame()
  write.csv(cor_mat, paste0(output_folder,"/VarCorrelations/02_CorMat/OverlapVME_",names(all_layers_overlap)[layer],"_cor.csv"))
  cor_long <- as.data.frame(cor_mat) %>%
    rownames_to_column(var = "var1") %>%
    pivot_longer(-var1, names_to = "var2", values_to = "cor")
  return(cor_long)
}) %>%
  set_names(names(all_layers_df)) %>%
  bind_rows(.id = "PeriodSSP") %>%
  mutate(
    period = ifelse(is.na(str_extract(PeriodSSP, "P\\d")), "baseline", str_extract(PeriodSSP, "P\\d")),
    ssp = ifelse(is.na(str_extract(PeriodSSP, "\\d-\\d\\.\\d")), "", str_extract(PeriodSSP, "\\d-\\d\\.\\d"))
  )

# Output correlation plots for overlapping cells ----
lapply(names(all_layers_overlap), function(x) {
  df <- filter(all_layers_overlap_cor, PeriodSSP == x)
  x_name <- ifelse(x == "baseline", "Reference", paste(df$period[1], "SSP", df$ssp[1]))
  p <- ggplot(data = df, aes(x = var1, y = var2, fill = cor)) +
    geom_tile(colour = "black") +
    geom_text(aes(label = ifelse(cor > 0.7 | cor < -0.7, "*", ""))) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", limit = c(-1,1), midpoint = 0) +
    scale_x_discrete(expand = c(0,0)) +
    scale_y_discrete(expand = c(0,0)) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text = element_text(size = 8),
          axis.line = element_blank(),
          axis.title = element_blank()) +
    labs(fill = "Correlation",
        title = x_name)
  ggsave(
    filename = paste0(output_folder,"/VarCorrelations/03_CorPlots/OverlapVME_",x,".jpg"),
    plot = p, width = 8, height = 6
  )
})

# Calculate correlation differences between adjacent time periods and baseline over all SSPs for full NAFO area ----
all_layers_cor_adj <- lapply(
  list.files(
    path = paste0(output_folder,"/VarCorrelations/02_CorMat"),
    pattern = "FullArea",
    full.names = TRUE
  ),
  function(x) read.csv(x, row.names = 1) |> as.matrix()
) %>%
  set_names(
    str_extract(
      list.files(path = paste0(output_folder,"/VarCorrelations/02_CorMat"), pattern = "FullArea"),
      "_(.+?)_", 
      group = 1
  ))

all_layers_cor_adj_dif <- lapply(ssp_all, function(sspoi) {
  layers_to_compare <- all_layers_cor_adj[c(1, grep(sspoi,names(all_layers_cor_adj)))]
  diff_layers <- list(
    layers_to_compare[[2]] - layers_to_compare[[1]],
    layers_to_compare[[3]] - layers_to_compare[[2]],
    layers_to_compare[[4]] - layers_to_compare[[3]],
    layers_to_compare[[5]] - layers_to_compare[[4]],
    layers_to_compare[[3]] - layers_to_compare[[1]],
    layers_to_compare[[4]] - layers_to_compare[[1]],
    layers_to_compare[[5]] - layers_to_compare[[1]]
  )
  names(diff_layers) <- c(
    "P1 - baseline",
    "P2 - P1",
    "P3 - P2",
    "P4 - P3",
    "P2 - baseline",
    "P3 - baseline",
    "P4 - baseline"
  )
  return(diff_layers)
}) %>%
  set_names(ssp_all) %>%
  unlist(recursive = FALSE)

all_layers_cor_adj_dif_df <- lapply(1:length(all_layers_cor_adj_dif), function(x) {
  write.csv(
    all_layers_cor_adj_dif[x],
    file = paste0(output_folder,"/VarCorrelations/04_CorDiffMat/FullArea_",names(all_layers_cor_adj_dif)[x],"_cordiff.csv")
  )
  cor_long <- as.data.frame(all_layers_cor_adj_dif[[x]]) %>%
    rownames_to_column(var = "var1") %>%
    pivot_longer(-var1, names_to = "var2", values_to = "cordiff")
  return(cor_long)
}) %>%
  set_names(names(all_layers_cor_adj_dif)) %>%
  bind_rows(.id = "DiffID")

# Output correlation difference plots ----
lapply(names(all_layers_cor_adj_dif), function(x) {
  df <- filter(all_layers_cor_adj_dif_df, DiffID == x)
  p <- ggplot(data = df, aes(x = var1, y = var2, fill = cordiff)) +
    geom_tile(colour = "black") +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
      limit = c(min(all_layers_cor_adj_dif_df$cordiff), max(all_layers_cor_adj_dif_df$cordiff)), midpoint = 0) +
    scale_x_discrete(expand = c(0,0)) +
    scale_y_discrete(expand = c(0,0)) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text = element_text(size = 8),
          axis.line = element_blank(),
          axis.title = element_blank()) +
    labs(fill = "Correlation difference",
         title = x)
  ggsave(
    filename = paste0(output_folder,"/VarCorrelations/05_CorDiffPlots/FullArea_",x,".jpg"),
    plot = p, width = 8, height = 6
  )
})

# Calculate correlation differences between adjacent time periods and baseline over all SSPs for overlapping PA cells ----
all_layers_cor_adj <- lapply(
  list.files(
    path = paste0(output_folder,"/VarCorrelations/02_CorMat"),
    pattern = "OverlapVME",
    full.names = TRUE
  ),
  function(x) read.csv(x, row.names = 1) |> as.matrix()
) %>%
  set_names(
    str_extract(
      list.files(path = paste0(output_folder,"/VarCorrelations/02_CorMat"), pattern = "OverlapVME"),
      "_(.+?)_", 
      group = 1
  ))

all_layers_cor_adj_dif <- lapply(ssp_all, function(sspoi) {
  layers_to_compare <- all_layers_cor_adj[c(1, grep(sspoi,names(all_layers_cor_adj)))]
  diff_layers <- list(
    layers_to_compare[[2]] - layers_to_compare[[1]],
    layers_to_compare[[3]] - layers_to_compare[[2]],
    layers_to_compare[[4]] - layers_to_compare[[3]],
    layers_to_compare[[5]] - layers_to_compare[[4]],
    layers_to_compare[[3]] - layers_to_compare[[1]],
    layers_to_compare[[4]] - layers_to_compare[[1]],
    layers_to_compare[[5]] - layers_to_compare[[1]]
  )
  names(diff_layers) <- c(
    "P1 - baseline",
    "P2 - P1",
    "P3 - P2",
    "P4 - P3",
    "P2 - baseline",
    "P3 - baseline",
    "P4 - baseline"
  )
  return(diff_layers)
}) %>%
  set_names(ssp_all) %>%
  unlist(recursive = FALSE)

all_layers_cor_adj_dif_df <- lapply(1:length(all_layers_cor_adj_dif), function(x) {
  write.csv(
    all_layers_cor_adj_dif[x],
    file = paste0(output_folder,"/VarCorrelations/04_CorDiffMat/OverlapVME_",names(all_layers_cor_adj_dif)[x],"_cordiff.csv")
  )
  cor_long <- as.data.frame(all_layers_cor_adj_dif[[x]]) %>%
    rownames_to_column(var = "var1") %>%
    pivot_longer(-var1, names_to = "var2", values_to = "cordiff")
  return(cor_long)
}) %>%
  set_names(names(all_layers_cor_adj_dif)) %>%
  bind_rows(.id = "DiffID")

# Output correlation difference plots ----
lapply(names(all_layers_cor_adj_dif), function(x) {
  df <- filter(all_layers_cor_adj_dif_df, DiffID == x)
  p <- ggplot(data = df, aes(x = var1, y = var2, fill = cordiff)) +
    geom_tile(colour = "black") +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
      limit = c(min(all_layers_cor_adj_dif_df$cordiff), max(all_layers_cor_adj_dif_df$cordiff)), midpoint = 0) +
    scale_x_discrete(expand = c(0,0)) +
    scale_y_discrete(expand = c(0,0)) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text = element_text(size = 8),
          axis.line = element_blank(),
          axis.title = element_blank()) +
    labs(fill = "Correlation difference",
         title = x)
  ggsave(
    filename = paste0(output_folder,"/VarCorrelations/05_CorDiffPlots/OverlapVME_",x,".jpg"),
    plot = p, width = 8, height = 6
  )
})
