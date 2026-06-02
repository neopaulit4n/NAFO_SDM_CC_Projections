
# Modelling with new methodology

# Load data ----
source("code/01_LoadData.R")

# Average predictors for each period and SSP
cmip_df_period_ssp <- cmip_df %>%
  filter(!is.na(period)) %>%
  group_by(lon, lat, period, ssp) %>%
  # summarise(across(sobavg:mldavg, \(x) mean(x, na.rm = TRUE)), .groups = "drop")  
  summarise(across(sobavg:mldavg, list(mean = ~mean(.x, na.rm = TRUE),
                                       min = ~min(.x, na.rm = TRUE),
                                       max = ~max(.x, na.rm = TRUE)),
                   .names = "{col}_{fn}")) %>%
  ungroup() %>%
  # Only keep max mldavg
  select(-(starts_with("mldavg") & !ends_with("_max")))


# Get study area extent and create spatial mask ----
sa <- terra::rast("data/raw/Bathy_Layers/GEBCO2024_FS005.tif") %>%
  terra::as.polygons(.) %>%
  sf::st_as_sf(.) %>%
  sf::st_transform(4326) %>%
  mutate(NRA_BNAM_b_tmp_mean = 1) %>%
  group_by(NRA_BNAM_b_tmp_mean) %>%
  summarise(geometry = sf::st_union(geometry))

sa_lims <- sf::st_coordinates(sa) %>% as.data.frame %>%
  select(X,Y)
sa_lims <- c(min(sa_lims$X), max(sa_lims$X),min(sa_lims$Y),max(sa_lims$Y))

# Transform CMIP data to raster
transform_cmip_to_raster <- function(data, varstat, period, ssp) {
  
  # Data must be ensembled dataframe with period and ssp columns
  
  # Select a layer of data to raster
  df <- data %>%
    filter(ssp == ssp, period == period) %>%
    select(lon, lat, !!sym(varstat)) %>%
    sf::st_as_sf(coords = c("lon","lat"), crs = 4326)
  
  # Transform into terra vector of points
  pts <- terra::vect(df)
  
  # Determine resolution of data in degrees
  # res <- round(cmip_ens_proj_df_period_cell$lon[2]-cmip_ens_proj_df_period_cell$lon[1], 5)
  res <- round(sort(unique(data$lon))[2] - sort(unique(data$lon))[1], 5)
  
  # Create template raster
  rast_template <- terra::rast(
    xmin = sa_lims[1], xmax = sa_lims[2], ymin = sa_lims[3], ymax = sa_lims[4],
    resolution = res,
    crs = "EPSG:4326"
  )
  # rast_template <- bathy_layers[[1]]
  
  # Rasterise points to grid
  rast_result <- terra::rasterize(pts, rast_template, field = varstat)
  
  # Resample to match BNAM raster
  rast_resamp <- terra::resample(rast_result, bnam_layers[[1]], method = "bilinear")
  return(rast_resamp)
  
}

vars <- colnames(select(cmip_df_period_ssp, contains("_")))

cmip_layers <- lapply(vars, transform_cmip_to_raster, 
                      data = cmip_df_period_ssp, 
                      period = "P1", 
                      ssp = "2-4.5")
names(cmip_layers) <- vars
# terra::plot(cmip_layers[[1]])

cmip_pred_df <- lapply(c(bathy_layers, cmip_layers), function(layer) {
  terra::extract(layer, 
                 select(resp_df, Start_Long_DD, Start_Lat_DD)) %>%
    select(-ID)
}) %>%
  bind_cols() %>%
  set_names(c(names(bathy_layers), names(cmip_layers)))


# Combine predictor and response dataframes ----
cmip_comb_df <- bind_cols(resp_df, cmip_pred_df) %>%
  mutate(VME_P_A = factor(VME_P_A, levels = c(0, 1), labels = c("Absence", "Presence"))) %>%
  drop_na()


# Getting variables ready ----
vme_group <- "black_corals"
vme_terrain_vars <- filter(terrain_topvars, VME_Group == vme_group) %>%
  pull(variable)
vme_vars <- c(vme_terrain_vars, names(cmip_layers))
vme_df <- filter(cmip_comb_df, VME_Group == vme_group) %>%
  select(all_of(c("VME_P_A", vme_vars)))


# Variable selection based on correlation and VIF ----
vars <- colnames(vme_df[,-1])

## Calculate correlation ----
cor_df <- cor(vme_df[, vars, drop = FALSE], 
                   method = "spearman",
                   use = "pairwise.complete.obs") %>%
  as.data.frame %>%
  pivot_longer(everything(), names_to = "var2", values_to = "cor") %>%
  mutate(var1 = rep(vars, each = length(vars))) %>%
  select(var1, var2, cor) #%>%
  # Assign numeric priorities according to Knudby et al 2013
  # mutate()

## Assign numeric priorities to variables ----
var_priority <- sapply(vars, function(v) {
  priority <- 0
  if (v %in% vme_terrain_vars) priority <- priority + 10000
  if (grepl("bstress|sobavg|tobavg|wobavg", v)) priority <- priority + 1000
  if (grepl("max|min", v)) priority <- priority + 200
  if (grepl("range", v)) priority <- priority + 100
  return(priority)
})
var_priority <- sort(var_priority, decreasing = TRUE)






removed_vars <- c()
equal_prio_vars <- c()
verbose <- TRUE
cor_matrix <- cor(vme_df[, vars, drop = FALSE], 
                  method = "spearman",
                  use = "pairwise.complete.obs")
for (i in 1:length(vars)) {
  var_i <- vars[i]
  if (verbose) cat("\nEvaluating variable:", var_i, "\n")
  
  # Skip if this variable was already removed
  if (var_i %in% removed_vars) {
    if (verbose) cat("  Skipped (already removed)\n")
    next
  }
  
  for (j in i:length(vars)) {
    var_j <- vars[j]
    if (verbose) cat("  Comparing with variable:", var_j, "\n")
    
    # Skip if same variable or already removed
    if (i == j || var_j %in% removed_vars) {
      if (verbose) cat("    Skipped (same variable or already removed)\n")
      next
    }
    
    # Check correlation/association
    if (abs(cor_matrix[var_i, var_j]) > 0.7) {
      
      # Remove the lower priority variable
      if (var_priority[var_i] > var_priority[var_j]) {
        removed_vars <- c(removed_vars, var_j)
        if (verbose) cat("    Removed", var_j, "due to high association with", var_i,
                         "(", round(cor_matrix[var_i, var_j], 2), ")\n")
      } else if (var_priority[var_i] == var_priority[var_j]) {
        equal_prio_vars <- rbind(equal_prio_vars, data.frame(var1 = var_i, var2 = var_j, cor = cor_matrix[var_i, var_j]))
        if (verbose) cat("    Equal priority for", var_i, "and", var_j,
                         "(", round(cor_matrix[var_i, var_j], 2), ")\n")
      } else {
        removed_vars <- c(removed_vars, var_i)
        if (verbose) cat("    Removed", var_i, "due to high association with", var_j,
                         "(", round(cor_matrix[var_i, var_j], 2), ")\n")
        break # No need to compare var_i with other variables if it's removed
      }
      
    }
  }
}
remaining_vars <- setdiff(vars, removed_vars)

ggplot(data = cor_df, aes(x = var1, y = var2, fill = cor)) +
  geom_tile(colour = "black") +
  geom_text(aes(label = ifelse(cor > 0.7 | cor < -0.7, "*", "")), size = 3) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", limit = c(-1,1), midpoint = 0) +
  scale_x_discrete(expand = c(0,0)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.line = element_blank()) +
  labs(title = paste("Correlation Matrix for VME Group:", vme_group), fill = "Correlation")

cor_df_equal <- filter(cor_df,
                       var1 %in% equal_prio_vars$var1 & var2 %in% equal_prio_vars$var2)
ggplot(data = cor_df_equal, aes(x = var1, y = var2, fill = cor)) +
  geom_tile(colour = "black") +
  geom_text(aes(label = ifelse(cor > 0.7 | cor < -0.7, "*", "")), size = 3) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", limit = c(-1,1), midpoint = 0) +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.line = element_blank()) +
  labs(title = paste("Correlation Matrix for VME Group:", vme_group), fill = "Correlation")


## Setting up GUI to eliminate remaining equal variables interactively ----

### Setup filenames
tmpfilename <- "temporary.txt"
outlistfilename <- "elimination_list.txt"
remainingfilename <- "remaining_list.csv"

vars_manually_removed <- c()
cor_df_equal_filtered <- filter(cor_df_equal, 
                                abs(cor) > 0.7,
                                !(var1 == var2)) %>%
  arrange(desc(abs(cor)))
row_comp <- 1

### Functions for buttons
write1 <- function() {
  tmpfile <- file(tmpfilename)
  writeLines(cor_df_equal_filtered$var1[1], tmpfile)
  close(tmpfile)
}
write2 <- function() {
  tmpfile <- file(tmpfilename)
  writeLines(cor_df_equal_filtered$var2[1], tmpfile)
  close(tmpfile)
}

### GUI widget
root <- tcltk::tktoplevel()
btn1 <- tcltk2::tk2button(root, text=paste(cor_df_equal_filtered$var1[1]), 
                          command = write1)
tcltk::tkpack(btn1)
btn2 <- tcltk2::tk2button(root, text=paste(cor_df_equal_filtered$var2[1]), 
                          command = write2)
tcltk::tkpack(btn2)
# Show correlation between the two variables in widget
tcltk::tkpack(tcltk::tklabel(root, text=paste("Correlation:", round(cor_df_equal_filtered$cor[1], 2))))

# Add comment section for user to fill in reasoning
tcltk::tkpack(tcltk::tklabel(root, text="Reason for elimination:"))
reason_entry <- tcltk::tkentry(root)

while (nrow(cor_df_equal_filtered)) { 
  
  # Name buttons after the two variables
  tcltk::tkconfigure(btn1, text=paste(cor_df_equal_filtered$var1[1]), command=write1)
  tcltk::tkconfigure(btn2, text=paste(cor_df_equal_filtered$var2[1]), command=write2)
  # Show correlation between the two variables in widget
  tcltk::tkconfigure(tcltk::tklabel(root, text=paste("Correlation:", round(cor_df_equal_filtered$cor[1], 2))))
  
  # Wait until the results have been written to the temporary file
  while (!file.exists(tmpfilename)) {a=1} # a is a dummy variable
  
  # Read the variable to be eliminated from the temporary file
  var2remove <- readLines(tmpfilename, 1) # Read one line
  file.remove(tmpfilename) # Delete the file
  
  # Eliminate variable by adding it to the list and removing it from 'correlations'
  vars_removed_manually <- c(vars_removed_manually, var2remove)
  cor_df_equal_filtered <- cor_df_equal_filtered %>%
    filter(!(var1 == var2remove | var2 == var2remove))
  
}

# Remove the GUI widget
tcltk::tkdestroy(root)

readLines("temporary.txt")





# Old code ----

pred_var_df <- vme_df[,-1]

select_vme_vars <- function(pred_var_df, vme_group, cor_threshold = 0.7, verbose = FALSE) {
  
  cat("\n==============================\n")
  cat("\nSelecting variables for VME group:", vme_group, "\n")
  
  vif_threshold_met <- FALSE
  cor_threshold_new <- cor_threshold
  
  while (!vif_threshold_met) {
    
    # Calculate Spearman correlations
    remove_cor_variables <- function(data, priority_list, current_threshold = 0.7, verbose = FALSE) {
      
      vars <- colnames(data)
      
      # Calculate correlation
      cor_matrix <- matrix(0, nrow = ncol(data), ncol = ncol(data))
      rownames(cor_matrix) <- colnames(cor_matrix) <- colnames(data)
      
      numeric_cor <- cor(data[, vars, drop = FALSE], 
                         method = "spearman",
                         use = "pairwise.complete.obs")
      cor_matrix[vars, vars] <- numeric_cor
      
      # Prioritise variables
      # priority_dict <- setNames(seq_along(priority_list), priority_list)
      
      # Knudby et al 2013 prioritisation
      # Prioritise max/min over range over mean
      # Prioritise max or min depending on variable
      # Prioritise annual over seasonal -> only for chlorophyll
      # Prioritise bottom over surface variables
      
      all_vars <- colnames(data)
      
      # Assign priorities (variables in priority list come first, by their order)
      var_priority <- sapply(all_vars, function(v) {
        if (v %in% priority_list) {
          return(priority_dict[v])
        } else {
          return(length(priority_list) + 1) # Lower priority for non-listed vars
        }
      })
      
      # Sort variables by priority
      sorted_vars <- all_vars[order(var_priority)]
      
      # Identify highly correlated pairs and remove lower priority variables
      removed_vars <- c()
      
      for (i in 1:length(sorted_vars)) {
        var_i <- sorted_vars[i]
        if (verbose) cat("\nEvaluating variable:", var_i, "\n")
        
        # Skip if this variable was already removed
        if (var_i %in% removed_vars) {
          if (verbose) cat("  Skipped (already removed)\n")
          next
        }
        
        for (j in i:length(sorted_vars)) {
          var_j <- sorted_vars[j]
          if (verbose) cat("  Comparing with variable:", var_j, "\n")
          
          # Skip if same variable or already removed
          if (i == j || var_j %in% removed_vars) {
            if (verbose) cat("    Skipped (same variable or already removed)\n")
            next
          }
          
          # Check correlation/association
          if (abs(cor_matrix[var_i, var_j]) > current_threshold) {
            # Remove the lower priority variable
            removed_vars <- c(removed_vars, var_j)
            if (verbose) cat("    Removed", var_j, "due to high association with", var_i,
                             "(", round(cor_matrix[var_i, var_j], 2), ")\n")
          }
        }
      }
      
      # Return remaining variables
      return(setdiff(colnames(data), removed_vars))
    }
    
    cat("\n   Removing correlated variables (threshold:", cor_threshold_new, ")\n")
    
    uncor_vars <- remove_cor_variables(pred_var_df, 
                                       priority_list = rf_prelim_imp %>%
                                         # filter(VME_Group == vme_group) %>%
                                         pull(Variable),
                                       current_threshold = cor_threshold_new,
                                       verbose = verbose)
    
    # Calculate VIF values
    cat("\n   Calculating VIF values for selected variables\n")
    
    corvif <- function(data, verbose = TRUE) {
      data <- as.data.frame(data)
      
      form    <- formula(paste("fooy ~ ", paste(strsplit(names(data), " "), collapse = " + ")))
      data   <- data.frame(fooy = 1 + rnorm(nrow(data)), data)
      lm_mod  <- lm(form, data)
      
      if (verbose) print(data.frame(vif=car::vif(lm_mod)))
      return(data.frame(vif=car::vif(lm_mod)))
    }
    vif_values <- corvif(pred_var_df[, uncor_vars], verbose = verbose)
    
    # Check if all VIF values are < 10
    if (all(vif_values$vif < 10)) {
      cat("\n     All VIF values are below threshold for VME group:", vme_group, "\n")
      vif_threshold_met <- TRUE
      # uncor_vars[[vme_group]] <- selected_vars_new
      # vars_vif[[vme_group]] <- vif_values_new
      return(list(selected_vars = uncor_vars,
                  vif_values = vif_values,
                  final_cor_threshold = cor_threshold_new))
    }
    else {
      cat("\n     VIF values exceed threshold (10); lowering correlation threshold and re-evaluating.\n")
      cor_threshold_new <- cor_threshold_new - 0.05
    }
  }
}

# selected_vme_vars <- lapply(unique(comb_df_compl$VME_Group), select_vme_vars, cor_threshold = 0.7) %>%
#   set_names(unique(comb_df_compl$VME_Group))

selected_vme_vars <- select_vme_vars(vme_df[,-1], vme_group = vme_group, cor_threshold = 0.7, verbose = TRUE)

# Create final dataframe with selected variables
vme_df_sel <- vme_df %>%
  select(all_of(c("VME_P_A", selected_vme_vars$selected_vars)))


