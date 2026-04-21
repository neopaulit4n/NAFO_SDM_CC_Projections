
# Exploring how variable correlations change across time periods (stationarity of variable correlations)

# Ellen suggested looking into Augmented Dickey-Fuller (ADF) or Wald tests but not sure if they're suitable for testing variable correlations changes?
# Examples I see for ADF are for only one variable time series, not correlations between variables

# I want to plot correlations between variable pairs over time

# baseline_layers <- as.list(vme_layers_baseline) %>%
#   # set_names(paste0("baseline.",names(vme_layers_baseline)))
#   set_names(names(vme_layers_baseline))
# baseline_layers <- c(baseline_layers)

baseline_layers <- c(cmip_layers, CHNETBL3 = bathy_layers$CHNETBL3) |>
  terra::rast()

# projection_layers <- unlist(cmip_layers_future, recursive = FALSE)
# projection_layers <- projection_layers[grep(paste(selected_vme_vars, collapse = "|"), names(projection_layers))]

# all_layers <- c(baseline_layers, projection_layers)

# projection_layers <- unlist(vme_layers_proj, recursive = FALSE)
projection_layers <- unlist(cmip_layers_proj, recursive = FALSE) |>
  lapply(X = _, 
    function(layer) {
      c(layer, CHNETBL3 = bathy_layers$CHNETBL3) |>
        terra::rast()
    })

all_layers <- list()
all_layers[[1]] <- baseline_layers
names(all_layers) <- "baseline"
all_layers <- c(all_layers, projection_layers)

all_layers_df <- lapply(all_layers, terra::as.data.frame)

# Output base dataframes used to calculate correlations between variables
if (!dir.exists("output/01_Exploratory/04_CorrelationMatrices")) dir.create("output/01_Exploratory/04_CorrelationMatrices")

lapply(1:length(all_layers_df), 
  function(layer) {
    write_csv(df, paste0("output/01_Exploratory/04_CorrelationMatrices/dataframes/",names(all_layers_df)[layer],"_df.csv"))
  })

# Let's extract the points where VME P/A data is - actually we don't want to do that for this
# vme_pred_df <- lapply(all_layers, function(layer) {
#   terra::extract(
#     layer, 
#     select(
#       filter(resp_df, VME_Group == vmeoi) %>%
#         drop_na(),
#       Start_Long_DD, Start_Lat_DD)
#     ) %>%
#     select(-ID)
# }) #%>%
#   bind_cols() %>%
#   set_names(names(all_layers))


# vme_pred_df <- terra::extract(vme_layers_current, select(resp_df, Start_Long_DD, Start_Lat_DD)) %>%
#   select(-ID)


# Calculate correlations
# vme_layers_cor <- cor(vme_pred_df, method = "spearman", use = "pairwise.complete.obs")
all_layers_cor <- lapply(1:length(all_layers_df), function(layer) {
  # cat(names(layer))
  cor_mat <- cor(all_layers_df[[layer]], method = "spearman", use = "complete.obs") |>
    as.data.frame()
  cor_long <- as.data.frame(cor_mat) %>%
    rownames_to_column(var = "var1") %>%
    pivot_longer(-var1, names_to = "var2", values_to = "cor")
  return(cor_long)
  # write.csv(cor_mat, paste0("output/01_Exploratory/04_CorrelationMatrices/correlation_matrices/",names(all_layers_df)[layer],"_cor.csv"))
}) %>%
  set_names(names(all_layers_df)) %>%
  bind_rows(.id = "PeriodSSP") %>%
  mutate(
    period = ifelse(is.na(str_extract(PeriodSSP, "P\\d")), "baseline", str_extract(PeriodSSP, "P\\d")),
    ssp = ifelse(is.na(str_extract(PeriodSSP, "\\d-\\d\\.\\d")), "", str_extract(PeriodSSP, "\\d-\\d\\.\\d"))
  )

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
    filename = paste0(output_folder,"/",vmeoi,"_plot_cor_AllCMIPVars_FullArea_",x_name,".jpg"),
    plot = p, width = 8, height = 6
  )
})


# Tidy dataframe
# vme_layers_cor_long <- vme_layers_cor %>%
#   as.data.frame(.) %>%
#   rownames_to_column(., var = "var1") %>%
#   pivot_longer(-var1, names_to = "var2", values_to = "cor") #%>%
#   # filter(var1 != var2)

# p1 <- ggplot(data = filter(vme_layers_cor, period == "baseline"), aes(x = var1, y = var2, fill = cor)) +
#   geom_tile(colour = "black") +
#   # geom_text(aes(label = ifelse(cor > 0.7 | cor < -0.7, "*", ""))) +
#   scale_fill_gradient2(low = "blue", mid = "white", high = "red", limit = c(-1,1), midpoint = 0) +
#   scale_x_discrete(expand = c(0,0)) +
#   scale_y_discrete(expand = c(0,0)) +
#   theme_classic() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1),
#         #axis.text = element_text(size = 6),
#         axis.line = element_blank(),
#         axis.title = element_blank()) +
#   labs(fill = "Correlation",
#        title = "Baseline")
# ggsave(paste0(main_output_folder,"cor_plot_baseline.jpg"), plot = p1)

# p2 <- ggplot(data = filter(vme_layers_cor, period != "baseline"), aes(x = var1, y = var2, fill = cor)) +
#   facet_wrap(~ssp + period) +
#   geom_tile(colour = "black") +
#   # geom_text(aes(label = ifelse(cor > 0.7 | cor < -0.7, "*", ""))) +
#   scale_fill_gradient2(low = "blue", mid = "white", high = "red", limit = c(-1,1), midpoint = 0) +
#   scale_x_discrete(expand = c(0,0)) +
#   scale_y_discrete(expand = c(0,0)) +
#   theme_classic() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1),
#         #axis.text = element_text(size = 6),
#         axis.line = element_blank(),
#         axis.title = element_blank()) +
#   labs(fill = "Correlation")
# ggsave(paste0(main_output_folder,"cor_plot_periodssp.jpg"), plot = p2)

# Testing equality of correlation matrices

# Use tool eqCorrMatTest from ldstatsHD packages (cran/ldstatsHD)
# ldstatsHD relies on camel dependency (cran/camel)

# ldstatsHD::eqCorrMatTest()


# EX2 <- ldstatsHD::pcorSimulatorJoint(nobs = 50, nclusters = 3, nnodesxcluster = c(40,40,40), 
#                           pattern = "pow", diffType = "cluster", dataDepend = "diag", 
#                           pdiff=0.5)	

# d1 <- as.matrix(vme_pred_df$`P1.1-2.6`)
# d2 <- as.matrix(vme_pred_df$`P1.2-4.5`)

# test1 <- ldstatsHD::eqCorrMatTest(
#   d1, 
#   d2, 
#   testStatistic = c("AS","max","exc"), 
#   excAdj = FALSE,
#   # testNullDist = "asyDep",
#   paired = TRUE, 
#   nite = 400)
