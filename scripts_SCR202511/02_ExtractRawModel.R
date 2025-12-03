
# Extract from raw 22 model data in chunks

# Extract all data and apply spatial mask to all dimensions
extract_data_constrained <- function(variable) {
  data <- ncvar_get(nc_file, variable, start = c(1,1,1,60,1), count = c(-1,-1,1,12,-1))
  data[data == -9999] <- NA
  
  # Apply spatial mask across all dimensions
  for (time_idx in 1:dim(data)[3]) {
    for (mod_idx in 1:dim(data)[4]) {
      slice <- data[,,time_idx,mod_idx]
      slice[!sa_mask] <- NA
      data[,,time_idx,mod_idx] <- slice
    }
  }
  return(data)
}

tosavg_data <- extract_data_constrained("tosavg")


