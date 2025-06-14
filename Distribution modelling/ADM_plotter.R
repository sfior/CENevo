### This script plots CEN present, past and future projections. Results are generated and imported from distributionModels.R

## Load packages
library("ggplot2")
library("gridExtra")
library("raster")
library("rasterVis")
library("RColorBrewer")

## Set parameters
plot_originals <- FALSE
cast <- "forecast"
#cast <- "hindcast"
if (cast == "hindcast") {
  cast.var.name <- "LGM"
} else if (cast == "forecast") {
  CMIP <- 6 # select CMIP 5 or 6
  if (CMIP==5) {
    cast.var.name <- "RCP45.2061-2080"
  } else if (CMIP==6) {
    cast.var.name <- "SSP370.2071-2100"
  }
}

## Define plotting extent
e_cen_plot <- extent(5, 12, 43.5, 48)

## Define the CRS (coordinate reference system). We define both longitude-latitude and a UTM (Zone 32) CRS.
prj_longlat <- "+init=epsg:4326"
prj_utm32 <- "+init=epsg:32632"

## Import results (output of Distribution_modeller.R)
setwd("/Users/hl636/Documents/Hirzi/ETHPHD/SDM/Dsylvestris/Data/")
if(cast == "forecast") {
  if (CMIP == 5) {
    pred.prob.current.masked.CentralAlps <- raster(paste0("Final_CEN_SDMs_forSimoneMS/", cast, "_CMIP", CMIP, "/pred.prob.current.masked.CentralAlps.2025-05-23_14_11_06.grd"))
    pred.prob.LGM.masked <- raster(paste0("Final_CEN_SDMs_forSimoneMS/", cast, "_CMIP", CMIP, "/pred.prob.LGM.masked.forecast.2025-05-23_14_11_06.grd"))
    pred.prob.LGM.masked.CentralAlps <- raster(paste0("Final_CEN_SDMs_forSimoneMS/", cast, "_CMIP", CMIP, "/pred.prob.LGM.masked.CentralAlps.forecast.2025-05-23_14_11_06.grd"))
  } else if (CMIP == 6) {
    pred.prob.current.masked.CentralAlps <- raster(paste0("Final_CEN_SDMs_forSimoneMS/", cast, "_CMIP", CMIP, "/pred.prob.current.masked.CentralAlps.2025-05-24_21_22_04.grd"))
    pred.prob.LGM.masked <- raster(paste0("Final_CEN_SDMs_forSimoneMS/", cast, "_CMIP", CMIP, "/pred.prob.LGM.masked.forecast.2025-05-24_21_22_04.grd"))
    pred.prob.LGM.masked.CentralAlps <-  raster(paste0("Final_CEN_SDMs_forSimoneMS/", cast, "_CMIP", CMIP, "/pred.prob.LGM.masked.CentralAlps.forecast.2025-05-24_21_22_04.grd"))
  }
} else if (cast == "hindcast") {
  pred.prob.current.masked.CentralAlps <- raster(paste0("Final_CEN_SDMs_forSimoneMS/", cast, "/pred.prob.current.masked.CentralAlps.2023-08-03_17_15_01.grd"))
  pred.prob.LGM.masked <- raster(paste0("Final_CEN_SDMs_forSimoneMS/", cast, "/pred.prob.LGM.masked.hindcast.2023-08-03_17_15_01.grd"))
  pred.prob.LGM.masked.CentralAlps <- raster(paste0("Final_CEN_SDMs_forSimoneMS/", cast, "/pred.prob.LGM.masked.CentralAlps.hindcast.2023-08-03_17_15_01.grd"))
}

## Import DEM time series
GLACIER_LANDSEA_DEM_TIMESERIES <-stack("Chelsa_TraCE_Glacier_LandSea_DEM_TimeSeries.grd")
GLACIER_LANDSEA_DEM_TIMESERIES <- crop(GLACIER_LANDSEA_DEM_TIMESERIES, e_cen_plot)
ENV_RASTERS <-stack("./ENV_RASTERS_TODAY_EUROPE_CHELSA_LARGE.grd") 
waterBodies_raster <- ENV_RASTERS$PH_5cm
waterBodies_raster <- crop(waterBodies_raster, e_cen_plot)

## Plot original results
if(plot_originals == TRUE) {
  # Define plotting parameters
  zeroCol <-"grey90"
  myPalette<-brewer.pal(11,"RdYlBu")
  my.at=seq(0, 1, by=0.1)
  my.brks=seq(0.05, 1.05, by=0.1)
  #myLabels <- c("NA (species not present)", "0.9 - low allele", "0.8 - low allele", "0.7 - low allele", "0.6 - low allele", "0.5 - high/low allele", "0.6 - high allele", "0.7 - high allele", "0.8 - high allele", "0.9 - high allele", "1.0 - high allele")
  myLabels <- c("NA (species not present)", "0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9", "1.0")
  myColorkey <- list(at=my.brks, labels=list(at=my.at, labels=myLabels))
  myTheme <- rasterTheme(region = c(zeroCol, myPalette))
  # Plot CEN distribution under present central Alpine lineage's mask 
  levelplot(pred.prob.current.masked.CentralAlps, par.settings = myTheme, at=my.at, colorkey=myColorkey, maxpixels = 1e7, margin = FALSE, main = "Continuous Spatial Projection (Masked Ensemble Model - Today)")
  # Plot CEN hindcast/forecast under future/past central Alpine lineage's distribution mask
  levelplot(pred.prob.LGM.masked, par.settings = myTheme, at=my.at, colorkey=myColorkey, maxpixels = 1e7, margin = FALSE, main = paste0("Continuous Spatial Projection (Masked Ensemble Model - ", cast.var.name, ")"))
  # Plot CEN hindcast/forecast under present central Alpine lineage's distribution mask
  levelplot(pred.prob.LGM.masked.CentralAlps, par.settings = myTheme, at=my.at, colorkey=myColorkey, maxpixels = 1e7, margin = FALSE, main = paste0("Continuous Spatial Projection (Masked Ensemble Model - ", cast.var.name, ")"))
}

## Plot with DEM baselayer
CEN_raster_stack <- stack(pred.prob.current.masked.CentralAlps, pred.prob.LGM.masked, pred.prob.LGM.masked.CentralAlps)
CEN_dem_raster_stack <- stack()
for(i in 1:nlayers(CEN_raster_stack)) {
  CEN_raster.dem <- CEN_raster_stack[[i]]
  CEN_raster.dem[CEN_raster.dem == 0.00] <- NA
  CEN_raster.dem <- CEN_raster.dem + 5000
  if(cast == "forecast" || i == 1 ) {
    base_map_dem_init <- GLACIER_LANDSEA_DEM_TIMESERIES[[211]]
  } else if (cast == "hindcast" && i >= 1) {
    base_map_dem_init <- GLACIER_LANDSEA_DEM_TIMESERIES[[1]]
  }
  base_map_dem_init[base_map_dem_init == -32768] <- NA # assign NA value
  CEN_raster.dem <- cover(CEN_raster.dem, base_map_dem_init)
  if(cast == "forecast" || i == 1 ) {
    CEN_raster.dem <- overlay(CEN_raster.dem, waterBodies_raster, fun=function(x, y) ifelse(is.na(y), NA, x))
  }
  proj4string(CEN_raster.dem) = CRS(prj_longlat)
  # Append to stack
  CEN_dem_raster_stack <- addLayer(CEN_dem_raster_stack, CEN_raster.dem)
}

## Define plotting parameters. Fixed z-range, so global min/max.
min_DEM_zValue <- min(minValue(GLACIER_LANDSEA_DEM_TIMESERIES)[which(minValue(GLACIER_LANDSEA_DEM_TIMESERIES) >-32768)])
max_DEM_zValue <- max(maxValue(GLACIER_LANDSEA_DEM_TIMESERIES))
breakpoints_DEM <- seq(min_DEM_zValue, max_DEM_zValue, length.out = 100)
breakpoints_projection <- seq(5000.05, 5001.05, by=0.1)
breakpoints <- c(breakpoints_DEM, breakpoints_projection)
color.palette=colorRampPalette(c('gray95','gray50'))
custom_colors <- c(color.palette(100), brewer.pal(11, "RdYlBu"))
my.at=seq(5000, 5001, by=0.1)
my.brks=seq(5000.05, 5001.05, by=0.1)
myLabels <- c("NA (species not present)", "0.9 - low allele", "0.8 - low allele", "0.7 - low allele", "0.6 - low allele", "0.5 - high/low allele", "0.6 - high allele", "0.7 - high allele", "0.8 - high allele", "0.9 - high allele", "1.0 - high allele")
myColorkey <- list(at=my.brks, labels=list(at=my.at, labels=myLabels))
# Plot
levelplot(CEN_dem_raster_stack[[1]], col.regions = custom_colors, at=breakpoints,  maxpixels = 1e8, margin = FALSE)
levelplot(CEN_dem_raster_stack[[2]], col.regions = custom_colors, at=breakpoints,  maxpixels = 1e8, margin = FALSE)
levelplot(CEN_dem_raster_stack[[3]], col.regions = custom_colors, at=breakpoints,  maxpixels = 1e8, margin = FALSE)
