# ------------------------------------------------------------------------------
# @Organization : Benjamin-Chung Lab
# @Project      : Environmental Risk Factors for Vivax Malaria in Peru
# @File         : fig2-map-study-area.R
# @Description  : Generate Figure 2: map of Loreto, Peru showing 31 study
#                 communities with annual P. vivax incidence rates.
# ------------------------------------------------------------------------------
rm(list=ls())

source(paste0(here::here(), "/0-config.R"))

library(sp)
library(raster)
library(terra)
library(sf)
library(tidyverse)
library(viridis)
library(RColorBrewer)
library(scales)
library(ggplot2)
library(ggspatial)

# Load spatial data -------------------------------------------------------
# read in Peru shapefile
districts = c("Punchana", "Alto Nanay", "San Juan Bautista", "Iquitos", "Belen")
peru_M_region_shp <- st_read(paste0(box_path_flame_erf, "shapefiles/per_admbnda_adm3_ign_20200714.shp")) %>% 
  filter(ADM2_ES == "Maynas") %>% filter(ADM3_ES %in% districts)

# Read in forest raster
lu_raster <- terra::rast(paste0(box_path_ucsf_shared, "rasters/lu_raster_epsg.tif"))
lu_raster[!(lu_raster %in% c(1, 26))] <- NA # keep only forest (1) and water (26)

# Aggregate raster to 100m x 100m resolution
lu_raster_agg <- terra::aggregate(lu_raster, fact=c(100/30, 100/30), fun="modal")

# Check and reproject raster if necessary
if (crs(lu_raster_agg) != st_crs(peru_M_region_shp)) {
  print("Reprojecting raster to match shapefile CRS...")
  lu_raster_agg <- terra::project(lu_raster_agg, st_crs(peru_M_region_shp)$wkt)
}

# Convert SpatRaster to a data frame for ggplot
lu_df <- as.data.frame(lu_raster_agg, xy = TRUE)
names(lu_df)[3] <- "value"

# Remove NA values and create separate datasets for forest and water
lu_df <- lu_df[!is.na(lu_df$value), ]
forest_df <- lu_df %>% filter(value == 1)
water_df <- lu_df %>% filter(value == 26)

# read in merged data for incidence calc
merged_df <- readRDS(paste0(box_path_flame_erf, 
                            "non-lagged-analysis-data_ext.RDS")) %>%
  filter(population > 15)

inc_df = merged_df %>% 
  group_by(comm_id, year, week) %>% 
  mutate(weekly_inc = sum(n_cases)/mean(population)*10000) %>% 
  ungroup() %>% 
  group_by(comm_id) %>% 
  summarize(mean_weekly_inc = mean(weekly_inc))

# centroids
comm_type_df = read.csv(paste0(box_path_flame_erf, "all_district_centroids_comm_type.csv")) %>%
  mutate(comm_id = sprintf("%03d", comm_id)) %>%
  mutate(comm_type = ifelse(comm_type == "", "Riverine", "Highway"))

centroids = read.csv(paste0(box_path_flame_erf, "centroids_final.csv")) %>%
  mutate(comm_id = sprintf("%03d", comm_id)) %>%
  left_join(comm_type_df, by = "comm_id")
inc_df = left_join(inc_df, centroids, by = "comm_id")

# calculate study area (km2)
centroids_matrix <- centroids %>%
  mutate(long = as.numeric(long), lat = as.numeric(lat)) %>%
  dplyr::select(long, lat)

# Get matrix of coordinates
coords <- as.matrix(centroids_matrix)

# Compute convex hull indices
hull_indices <- chull(coords[, 1], coords[, 2])

# Use convex hull to order the points and close the polygon
hull_coords <- coords[c(hull_indices, hull_indices[1]), ]  # Close the polygon

# Create the convex hull polygon as an sf object
hull_polygon <- st_polygon(list(hull_coords)) |> st_sfc(crs = 4326)

# Calculate the area (in square kilometers)
area_km2 <- st_area(hull_polygon) / 1e6

inc_palette_reversed <- plasma(9, direction = -1)


ggplot() +
  geom_sf(data = peru_M_region_shp, fill = NA, color = "black") +
  ggspatial::annotation_map_tile(type = "osm", zoom = 10, alpha = 0.8)

# Create ggplot map
map_study_area <- ggplot() +
  
  # Add background map tiles
  ggspatial::annotation_map_tile(type = "osm", zoom = 10, alpha = 0.8) +
  
  # Start with shapefile to establish bounds
  geom_sf(data = peru_M_region_shp, fill = NA, color = "black") +
  
  # Add forest areas
  geom_raster(data = forest_df, aes(x = x, y = y), fill = "chartreuse4", alpha = 0.9) +
  
  # Add water areas  
  geom_raster(data = water_df, aes(x = x, y = y), fill = "dodgerblue4", alpha = 0.9) +
  
  # Add community points with incidence coloring and different shapes
  geom_point(data = inc_df,
             aes(x = long, y = lat, fill = mean_weekly_inc, shape = comm_type),
             color = "black",
             size = 3.5,
             stroke = 0.5,
             position = position_jitter(width = 0.01, height = 0.01, seed = 123)) +
  
  # Define custom shapes (squares for highway, circles for riverine)
  scale_shape_manual(values = c("Highway" = 22,    # filled square
                                "Riverine" = 21),   # filled circle
                     name = "Community Type") +
  
  # Color scale for incidence
  scale_fill_gradientn(colors = inc_palette_reversed,
                       name = "Malaria Incidence\nper 10,000 person-weeks",
                       na.value = "transparent") +
  
  # Set plot limits to raster extents
  coord_sf(xlim = c(-74.5, -73), ylim = c(-4.2, -3.0), expand = FALSE) +
  
  theme_void() +
  theme(
    legend.position = "bottom",
    legend.key.width = unit(1, "cm"),
    legend.key.height = unit(0.5, "cm"),
    legend.margin = margin(10, 10, 10, 10),
    legend.text = element_text(family = "serif", size = 10),
    legend.title = element_text(family = "serif", size = 10),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(10, 10, 10, 10)
  )

# Save the map
ggsave(paste0(figure_path, "figure-2.pdf"),
       plot = map_study_area,
       width = 8, height = 6, units = "in", dpi = 300)
