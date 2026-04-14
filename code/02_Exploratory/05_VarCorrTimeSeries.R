
# Exploring how variable correlations change across time periods (stationarity of variable correlations)

# Ellen suggested looking into Augmented Dickey-Fuller (ADF) or Wald tests but not sure if they're suitable for testing variable correlations changes?
# Examples I see for ADF are for only one variable time series, not correlations between variables

# I want to plot correlations between variable pairs over time


# Which layers should we calculate correlations for?
selected_vme_vars

baseline_layers <- 
baseline_layers <- as.list(vme_layers_current) %>%
  # set_names(paste0("baseline.",names(vme_layers_current)))
  set_names(names(vme_layers_current))
baseline_layers <- c(baseline_layers)

# projection_layers <- unlist(cmip_layers_future, recursive = FALSE)
# projection_layers <- projection_layers[grep(paste(selected_vme_vars, collapse = "|"), names(projection_layers))]

# all_layers <- c(baseline_layers, projection_layers)

projection_layers <- unlist(vme_layers_future, recursive = FALSE)

all_layers <- list()
all_layers[[1]] <- vme_layers_current
names(all_layers) <- "baseline"
all_layers <- c(all_layers, projection_layers)

# Let's extract the points where VME P/A data is
vme_pred_df <- lapply(all_layers, function(layer) {
  terra::extract(layer, select(resp_df, Start_Long_DD, Start_Lat_DD)) %>%
    select(-ID)
}) #%>%
#   bind_cols() %>%
#   set_names(names(all_layers))


# vme_pred_df <- terra::extract(vme_layers_current, select(resp_df, Start_Long_DD, Start_Lat_DD)) %>%
#   select(-ID)


# Calculate correlations
# vme_layers_cor <- cor(vme_pred_df, method = "spearman", use = "pairwise.complete.obs")
vme_layers_cor <- lapply(vme_pred_df, function(layer) {
  cor_mat <- cor(layer, method = "spearman", use = "pairwise.complete.obs")
  cor_long <- as.data.frame(cor_mat) %>%
    rownames_to_column(var = "var1") %>%
    pivot_longer(-var1, names_to = "var2", values_to = "cor")
  return(cor_long)
}) %>%
  bind_rows(.id = "PeriodSSP") %>%
  mutate(
    period = ifelse(is.na(str_extract(PeriodSSP, "P\\d")), "baseline", str_extract(PeriodSSP, "P\\d")),
    ssp = ifelse(is.na(str_extract(PeriodSSP, "\\d-\\d\\.\\d")), "", str_extract(PeriodSSP, "\\d-\\d\\.\\d"))
  )


# Tidy dataframe
# vme_layers_cor_long <- vme_layers_cor %>%
#   as.data.frame(.) %>%
#   rownames_to_column(., var = "var1") %>%
#   pivot_longer(-var1, names_to = "var2", values_to = "cor") #%>%
#   # filter(var1 != var2)

p1 <- ggplot(data = filter(vme_layers_cor, period == "baseline"), aes(x = var1, y = var2, fill = cor)) +
  geom_tile(colour = "black") +
  # geom_text(aes(label = ifelse(cor > 0.7 | cor < -0.7, "*", ""))) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", limit = c(-1,1), midpoint = 0) +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        #axis.text = element_text(size = 6),
        axis.line = element_blank(),
        axis.title = element_blank()) +
  labs(fill = "Correlation",
       title = "Baseline")
ggsave(paste0(main_output_folder,"cor_plot_baseline.jpg"), plot = p1)

p2 <- ggplot(data = filter(vme_layers_cor, period != "baseline"), aes(x = var1, y = var2, fill = cor)) +
  facet_wrap(~ssp + period) +
  geom_tile(colour = "black") +
  # geom_text(aes(label = ifelse(cor > 0.7 | cor < -0.7, "*", ""))) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", limit = c(-1,1), midpoint = 0) +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        #axis.text = element_text(size = 6),
        axis.line = element_blank(),
        axis.title = element_blank()) +
  labs(fill = "Correlation")
ggsave(paste0(main_output_folder,"cor_plot_periodssp.jpg"), plot = p2)
