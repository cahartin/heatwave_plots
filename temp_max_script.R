# Climate Data Contour Map Script
# This script reads NetCDF climate data, subsets to a region,
# calculates temporal average, and creates a contour map

# Load required libraries
library(ncdf4)
library(ggplot2)
library(reshape2)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)

# Input file path
# folder name
var <- "Tw" #Tw, Ta

# Variable name 
temp_var <- "Tw"  #Tw, tas

# Degrees
degree <- "4deg" #2deg, 3deg, 4deg

# Sats
stat <- "mean" #mean, 3hrmax
fun <- "mean" #mean, max

region <- "Bang" #UAE, Bang
#### MAKE SURE TO CHANGE REGION IN RASTER STEP 2 ###


# directory
directory <- "C:/Users/corin/OneDrive/Documents/GitHub/heatwave_plots"

nc_files <- list.files(path = paste0(directory, "/", degree,"/",var), 
                       pattern = "*.nc", full.names = TRUE)

# Region of interest - Bangladesh
lon_min <- 87  # Western boundary
lon_max <- 94   # Eastern boundary
lat_min <- 19    # Southern boundary
lat_max <- 27     # Northern boundary

Bang <- c(87, 94, 19, 27)

# Region of interest - UAE
#lon_min <- 45  # Western boundary
#lon_max <- 60   # Eastern boundary
#lat_min <- 20    # Southern boundary
#lat_max <- 30     # Northern boundary

UAE <- c(45, 60, 20, 30)

# ---------------------------------------------------
# STEP 1: Process each model file
# ---------------------------------------------------

cat("\n========================================\n")
cat("Processing", length(nc_files), "model files\n")
cat("========================================\n\n")

# Initialize lists to store data from each model
model_data_list <- list()
lon_subset <- NULL
lat_subset <- NULL

# Loop through each model file
for (i in seq_along(nc_files)) {
  nc_file <- nc_files[i]
  cat("Processing model", i, "of", length(nc_files), ":", basename(nc_file), "\n")
  
  # Check if file exists
  if (!file.exists(nc_file)) {
    warning("File not found: ", nc_file, " - skipping")
    next
  }
  
  # Open NetCDF file
  nc_data <- nc_open(nc_file)
  
  units <- ncatt_get(nc_data, temp_var, "units")
  
  print(units$value)
  
  # Read full coordinate dimensions first (these are small)
  lon <- ncvar_get(nc_data, "lon")  # or "longitude"
  lat <- ncvar_get(nc_data, "lat")  # or "latitude"
  
  if (i == 1) {
    cat("Full data range - Lon:", range(lon), "Lat:", range(lat), "\n")
  }
  
  # Find indices for the region of interest
  lon_idx <- which(lon >= lon_min & lon <= lon_max)
  lat_idx <- which(lat >= lat_min & lat <= lat_max)
  
  # Check if we found valid indices
  if (length(lon_idx) == 0 | length(lat_idx) == 0) {
    warning("No data found in specified region for this model - skipping")
    nc_close(nc_data)
    next
  }
  
  if (i == 1) {
    cat("Subsetting to", length(lon_idx), "longitude points and", 
        length(lat_idx), "latitude points\n")
    # Store coordinates from first model
    lon_subset <- lon[lon_idx]
    lat_subset <- lat[lat_idx]
  }
  
  # Get variable info to determine dimensions
  var_info <- nc_data$var[[temp_var]]
  ndims <- var_info$ndims
  
  # Determine the dimension order and prepare start/count vectors
  dim_names <- sapply(var_info$dim, function(x) x$name)
  
  if (i == 1) {
    cat("Variable dimensions:", dim_names, "\n")
  }
  
  # Initialize start and count vectors
  start_vec <- rep(1, ndims)
  count_vec <- rep(-1, ndims)  # -1 means read all
  
  # Find which dimensions are lon, lat, time
  lon_dim_idx <- which(dim_names %in% c("lon", "longitude", "x"))
  lat_dim_idx <- which(dim_names %in% c("lat", "latitude", "y"))
  time_dim_idx <- which(dim_names %in% c("time", "t"))
  
  # Set start and count for longitude
  if (length(lon_dim_idx) > 0) {
    start_vec[lon_dim_idx] <- min(lon_idx)
    count_vec[lon_dim_idx] <- length(lon_idx)
  }
  
  # Set start and count for latitude
  if (length(lat_dim_idx) > 0) {
    start_vec[lat_dim_idx] <- min(lat_idx)
    count_vec[lat_dim_idx] <- length(lat_idx)
  }
  
  # Read only the subsetted data in R
  temp_subset <- ncvar_get(nc_data, temp_var, start = start_vec, count = count_vec)
  
  cat("Subsetted data dimensions:", dim(temp_subset), "\n")
  
  # Close the NetCDF file
  nc_close(nc_data)
  
  # Calculate temporal max for this model
temp <- apply(temp_subset, c(1, 2), fun, na.rm = TRUE)

  # Store the temporal average for this model
  model_data_list[[i]] <- temp
  
  cat("Model", i, "processed successfully\n\n")
}

# Remove NULL entries (from skipped files)
model_data_list <- model_data_list[!sapply(model_data_list, is.null)]

# Check if we have any valid models
if (length(model_data_list) == 0) {
  stop("No valid model data found. Check your file paths and settings.")
}

cat("Successfully processed", length(model_data_list), "models\n")

# Calculate ensemble average across all models
cat("\nCalculating ensemble average across", length(model_data_list), "models\n")

# Stack all model data into a 3D array (lon x lat x models)
model_array <- abind::abind(model_data_list, along = 3)

# Calculate ensemble mean across the model dimension (dimension 3)
# # convert from K to C
constant_value <- 273.15
ensemble_avg <- apply(model_array, c(1, 2), mean, na.rm = TRUE) - constant_value

cat("final data dimensions:", dim(ensemble_avg), "\n")

cat("Ensemble average dimensions:", dim(ensemble_avg), "\n")
cat("Ensemble average range:", range(ensemble_avg, na.rm = TRUE), "\n")

# Also calculate ensemble standard deviation
ensemble_sd <- apply(model_array, c(1, 2), sd, na.rm = TRUE)

# ---------------------------------------------------
# STEP 2: save data as a raster
# ---------------------------------------------------
cat("\nsaving data as a raster...\n")

library(terra)
ensemble_avg_r <- rast(aperm(ensemble_avg))
ext(ensemble_avg_r) <- Bang
crs(ensemble_avg_r) <- "EPSG:4326"
writeRaster(ensemble_avg_r, file = paste0(region, "_", degree, "_", temp_var, "_", stat, "_", ".tiff"), overwrite=TRUE)

# ---------------------------------------------------
# STEP 3: test figure
# ---------------------------------------------------
cat("\ntest figure...\n")

plot(ensemble_avg_r)

cat("\nHave a nice day!\n")