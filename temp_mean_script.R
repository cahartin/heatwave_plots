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
var <- "Tw"

# Variable name 
temp_var <- "Tw"  

# Degrees
degree <- "2deg"

# Sats
stat <- "mean"

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

# Output files
output_file_mean <- paste0(degree, "_", temp_var, "_", stat, "_", "climate_contour_map_ensemble.png")
output_file_range <- paste0(degree, "_", temp_var, "_", stat, "_", "climate_contour_map_ensemble_range.png")


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
  
  # Calculate temporal mean for this model
  temp_mean <- apply(temp_subset, c(1, 2), mean, na.rm = TRUE) 
  
  # Store the temporal average for this model
  model_data_list[[i]] <- temp_mean
  
  cat("Model", i, "processed successfully\n\n")
}

# Remove NULL entries (from skipped files)
model_data_list <- model_data_list[!sapply(model_data_list, is.null)]

# Check if we have any valid models
if (length(model_data_list) == 0) {
  stop("No valid model data found. Check your file paths and settings.")
}

cat("Successfully processed", length(model_data_list), "models\n")

# ---------------------------------------------------
# STEP 2: Calculate multi-model ensemble statistics
# ---------------------------------------------------

cat("\n========================================\n")
cat("Calculating ensemble statistics\n")
cat("========================================\n\n")

# Convert list to 3D array [lon, lat, model]
model_array <- array(unlist(model_data_list), 
                     dim = c(dim(model_data_list[[1]]), length(model_data_list)))

# convert from K to C
constant_value <- 273.15

cat("Model array dimensions:", dim(model_array), "\n")

# Calculate ensemble mean (average across models)
ensemble_mean <- temp_mean
ensemble_mean <- ensemble_mean - constant_value
cat("Ensemble mean calculated\n")

# Calculate ensemble range (min and max)
ensemble_min <- apply(model_array, c(1, 2), min, na.rm = TRUE)
ensemble_min <- ensemble_min - constant_value
ensemble_max <- apply(model_array, c(1, 2), max, na.rm = TRUE)
ensemble_max <- ensemble_max - constant_value
ensemble_range <- ensemble_max - ensemble_min
cat("Ensemble range calculated\n")


# ---------------------------------------------------
# STEP 3: save data for plotting
# ---------------------------------------------------
cat("\nsaving data as a raster...\n")

library(terra)
ensemble_mean_r <- rast(aperm(ensemble_mean))
ext(ensemble_mean_r) <- Bang
crs(ensemble_mean_r) <- "EPSG:4326"
writeRaster(ensemble_mean_r, file = paste0(degree, "_", temp_var, "_", stat, "_", ".tiff"), overwrite=TRUE)

# # ---------------------------------------------------
# # STEP 4: Prepare data for plotting
# # ---------------------------------------------------
#
 cat("\nPreparing ensemble data for plotting...\n")

 # Create a data frame for ensemble mean
 grid_mean <- expand.grid(lon = lon_subset, lat = lat_subset)
 grid_mean$temp <- as.vector(ensemble_mean)
 grid_mean <- grid_mean[!is.na(grid_mean$temp), ]


 # Prepare range data
 grid_range <- expand.grid(lon = lon_subset, lat = lat_subset)
 grid_range$temp <- as.vector(ensemble_range)
 grid_range <- grid_range[!is.na(grid_range$temp), ]

 # ---------------------------------------------------
 # STEP 4: Create contour maps
 # ---------------------------------------------------

 # Get country boundaries (shared across all plots)
 world <- ne_countries(scale = "medium", returnclass = "sf")
 world_cropped <- st_crop(world, xmin = lon_min, xmax = lon_max,
                          ymin = lat_min, ymax = lat_max)

 ## PLOT 1: Ensemble Mean
 cat("Creating mean ensemble contour map...\n")

 p_mean <- ggplot() +
   geom_contour_filled(data = grid_mean, aes(x = lon, y = lat, z = temp), bins = 20) +
   geom_contour(data = grid_mean, aes(x = lon, y = lat, z = temp),
                color = "black", alpha = 0.3) +
   geom_sf(data = world_cropped, fill = NA, color = "black", linewidth = 0.5) +
   coord_sf(xlim = c(lon_min, lon_max), ylim = c(lat_min, lat_max), expand = FALSE) +
   labs(
     title = paste0(degree, "_", temp_var),
     subtitle = paste0("Region: ", lon_min, "°E to ", lon_max, "°E, ",
                       lat_min, "°N to ", lat_max, "°N"),
     x = "Longitude",
     y = "Latitude",
     fill = "WBGT\n (°C)"
   ) +
   theme_minimal() +
   theme(
     plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
     plot.subtitle = element_text(hjust = 0.5, size = 10),
    legend.position = "right",
     panel.background = element_rect(fill = "white"),
     panel.grid.major = element_line(color = "gray90", linewidth = 0.3)
   ) +
   scale_fill_viridis_d(option = "plasma")

# ggsave(output_file_mean, plot = p_mean, width = 10, height = 8, dpi = 300)
 cat("Ensemble mean plot saved to:", output_file_mean, "\n")
 print(p_mean)
#
# ## PLOT 2: Model Spread/Range
# cat("\nCreating model spread/range map...\n")
#
# p_range <- ggplot() +
#   geom_contour_filled(data = grid_range, aes(x = lon, y = lat, z = temp), bins = 15) +
#   geom_contour(data = grid_range, aes(x = lon, y = lat, z = temp),
#                color = "black", alpha = 0.3, bins = 15) +
#   geom_sf(data = world_cropped, fill = NA, color = "black", linewidth = 0.5) +
#   coord_sf(xlim = c(lon_min, lon_max), ylim = c(lat_min, lat_max), expand = FALSE) +
#   labs(
#     title = paste0(degree, "_", temp_var, "_", "Range (Max - Min)"),
#     subtitle = paste0("Region: ", lon_min, "°E to ", lon_max, "°E, ",
#                       lat_min, "°N to ", lat_max, "°N"),
#     x = "Longitude",
#     y = "Latitude",
#     fill = "Range\n(°C)"
#   ) +
#   theme_minimal() +
#   theme(
#     plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
#     plot.subtitle = element_text(hjust = 0.5, size = 10),
#     legend.position = "right",
#     panel.background = element_rect(fill = "white"),
#     panel.grid.major = element_line(color = "gray90", linewidth = 0.3)
#   ) +
#   scale_fill_viridis_d(option = "cividis")
#
# ggsave(output_file_range, plot = p_range, width = 10, height = 8, dpi = 300)
# cat("Model spread plot saved to:", output_file_range, "\n")
# print(p_range)


cat("\nHave a nice day!\n")
