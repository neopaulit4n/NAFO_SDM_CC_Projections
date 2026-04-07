# Map of habitat suitability using variable partial dependence data - troubleshooting why projections seem weird
# Weighting by variable importance?

cmip_all <- bind_rows(
  mutate(current_df, period = "Baseline", ssp = NA),
  cmip_df_period_ssp
) %>%
  select(-starts_with("month_"))

varstat <- "tosavg_range"
p <- lapply(unique(cmip_all$period), function(poi) {
  lapply(unique(cmip_all$ssp), function(sspoi) {
    data <- filter(cmip_all, period == poi, ssp %in% sspoi)
    data <- cmip_all
    ggplot() +
      theme_bw() +
      geom_tile(data = data, aes(x = lon, y = lat, fill = !!sym(varstat))) +
      facet_wrap(~ period + ssp) +
      cmocean::scale_fill_cmocean(varstat, name = "thermal",
        limits = c(min(data[[varstat]]), max(data[[varstat]]))) +
      labs(title = varstat,
        x = "Longitude",
        y = "Latitude") #+
      # theme(legend.position = "none")      
  })
})


