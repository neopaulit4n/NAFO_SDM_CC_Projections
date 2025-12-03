 
# stars method

library(tidyverse)
library(stars)

# Read and ensemble-average
nc_data <- read_ncdf("data/nafo_fishingfoot_2015_2099_cmip22_sorall_identity.nc",
                     var = "tosavg", proxy = TRUE)
nc_ensemble_mean <- nc_data %>%
  st_apply(MARGIN = c("lon", "lat", "time", "lev"), 
           FUN = mean, 
           na.rm = TRUE,
           .fname = "ensemble_mean")

# Save to new netCDF file
write_stars(nc_ensemble_mean, "data/nafo_fishingfoot_2015_2099_cmip22_sorall_identity_tosavg_mavg.nc")

# Later, load the smaller file
nc_avg <- read_stars("ensemble_averaged_data.nc")


# Read the netCDF file
nc_tosavg <- read_ncdf("data/nafo_fishingfoot_2015_2099_cmip22_sorall_identity.nc",
                     var = "tosavg",
                     # ncsub = cbind(start = c(1, 1, 1, 1, 1),
                     #               count = c(115, 84, 1, 1, 1)),
                     proxy = TRUE)

nc_tosavg_ens <- nc_tosavg %>%
  st_apply(MARGIN = c("lon", "lat", "time", "lev"), 
           FUN = mean, 
           na.rm = TRUE,
           .fname = "ensemble_mean")

# nc_tosavg <- filter(nc_tosavg, ssp == 126, ens == 1, time == 1008408)
nc_tosavg <- read_stars(nc_tosavg_ens, proxy = FALSE)

# Create time periods 
start_date <- as.Date("2020-01-15")  # adjust as needed
time_periods <- data.frame(
  time_index = 1:1020,
  date = seq(start_date, by = "month", length.out = 1020)
) %>%
  mutate(
    year = year(date),
    period = case_when(
      year >= 2020 & year <= 2039 ~ "P1",
      year >= 2040 & year <= 2059 ~ "P2",
      year >= 2060 & year <= 2079 ~ "P3",
      year >= 2080 & year <= 2099 ~ "P4"
    )
  )

# Add period information to your stars object
nc_data <- nc_data %>%
  mutate(period = time_periods$period[time])

# Calculate summary statistics by period
summary_stats <- nc_data %>%
  group_by(period, ssp, ens) %>%
  summarise(
    mean = mean(.data[[1]], na.rm = TRUE),
    min = min(.data[[1]], na.rm = TRUE),
    max = max(.data[[1]], na.rm = TRUE),
    sd = sd(.data[[1]], na.rm = TRUE),
    .groups = "keep"
  )



# Use ncdf4 to create an ensembled ncdf file
library(ncdf4)

# create_ensemble_averaged_netcdf <- function(input_file, output_file) {
  # Open original file
nc_orig <- nc_open("data/nafo_fishingfoot_2015_2099_cmip22_sorall_identity.nc")

# Get dimensions (excluding ensemble)
dims_to_keep <- c("lon", "lat", "time", "lev")  # adjust names as needed
new_dims <- list()

for(dim_name in dims_to_keep) {
  dim_vals <- ncvar_get(nc_orig, dim_name)
  dim_units <- ncatt_get(nc_orig, dim_name, "units")$value
  new_dims[[dim_name]] <- ncdim_def(dim_name, dim_units, dim_vals)
}

# Process variables
var_defs <- list()
var_data <- list()

for (var_name in names(nc_orig$var)) {
  # Get original variable info
  orig_var <- nc_orig$var[[var_name]]
  units <- ncatt_get(nc_orig, var_name, "units")$value
  
  # Read and average data
  data_5d <- ncvar_get(nc_orig, var_name)
  data_4d <- apply(data_5d, MARGIN = c(1,2,3,4), FUN = mean, na.rm = TRUE)
  
  # Create variable definition
  var_defs[[var_name]] <- ncvar_def(
    name = var_name,
    units = units,
    dim = new_dims,
    missval = orig_var$missval
  )
  
  var_data[[var_name]] <- data_4d
  cat("Processed variable:", var_name, "\n")
}

# Create new file
nc_new <- nc_create(output_file, var_defs)

# Write data
for(var_name in names(var_data)) {
  ncvar_put(nc_new, var_name, var_data[[var_name]])
}

# Copy relevant global attributes
global_atts <- ncatt_get(nc_orig, 0)
for(att_name in names(global_atts)) {
  ncatt_put(nc_new, 0, att_name, global_atts[[att_name]])
}

# Add processing note
ncatt_put(nc_new, 0, "processing_note", 
          "Ensemble mean calculated across all ensemble members")

nc_close(nc_orig)
nc_close(nc_new)

# cat("Ensemble-averaged file created:", output_file, "\n")
# }

# Use the function
create_ensemble_averaged_netcdf("your_file.nc", "ensemble_averaged_data.nc")



