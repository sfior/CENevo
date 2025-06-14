### This script plots and compares CEN ADM forecasts+shifts with gradient forest forecasts+shifts
# Run the following scripts first: gradientForest_glacialGenomicOffset.R; gradientForest_genomicOffset.R; ADM_plotter.R
plot_originals <- FALSE

# Define cropping extent and mask
e_mask <- extent(5, 12, 43.5, 48)
current.masked.cropped <- crop(current.mask, extent(e_mask))

# Calculate the difference between present-day and forecasted CEN projections.
deltaCEN <- pred.prob.LGM.masked.CentralAlps - pred.prob.current.masked.CentralAlps
deltaCEN <- crop(deltaCEN, extent(e_mask))
deltaCEN <- mask(deltaCEN, current.masked.cropped, maskvalue=0, updatevalue = NA)
deltaCEN <- -deltaCEN

# Calculate the multivariate Euclidean distance between present-day and forecast genome-wide GF projections.
genomicOffset.df <- as.data.frame((rowSums((CLIM.alps.future.df.GF[,-c(1,2,3,4)] - CLIM.alps.df.GF[,-c(1,2,3,4)])^2))^0.5)
# Make dataframe of x,y,dist
genomicOffset.df <- data.frame(CLIM.alps.df.GF[,c(1,2)], genomicOffset.df)
# Convert to raster
genomicOffset <- rasterFromXYZ(genomicOffset.df)
crs(genomicOffset) <- crs(CLIM.alps)
# Crop and mask projections by SDM layer
genomicOffset <- crop(genomicOffset, extent(e_mask))
genomicOffset <- mask(genomicOffset, current.masked.cropped, maskvalue=0, updatevalue = NA)

# Import the future flacial genomic offset 
glacialGenomicOffset <- crop(CLIM.alps.future.GF.dist.rast, extent(e_mask))
glacialGenomicOffset <- mask(glacialGenomicOffset, current.masked.cropped, maskvalue=0, updatevalue = NA)

# # Test correlation
# x <- deltaCEN
# #x <- genomicOffset
# y <- genomicOffset
# #y <- glacialGenomicOffset
# cor(x@data@values, y@data@values, method = "spearman", use = "pairwise.complete.obs") # this is the same as ENMTools' raster.cor function
# cor.test(x@data@values, y@data@values, method = "spearman", use = "pairwise.complete.obs")

## Plot results (simple plot; for colour legend)
if(plot_originals == TRUE) {
  # Define plotting parameters
  nbreaks <- 11
  zeroCol <-"grey90"
  myPalette <- rev(brewer.pal(nbreaks,"RdYlBu"))
  my.at.CEN <- seq(-0.2, 0.37, length.out = nbreaks)
  #my.at.GF <- seq(min(genomicOffset@data@min, glacialGenomicOffset@data@min), max(genomicOffset@data@max, glacialGenomicOffset@data@max), length.out = nbreaks)
  my.at.GF <- seq(0, 0.038, length.out = nbreaks)
  my.brks.CEN <- my.at.CEN + ((my.at.CEN[2] - my.at.CEN[1])/2)
  my.brks.GF <- my.at.GF + ((my.at.GF[2] - my.at.GF[1])/2)
  myLabels.CEN <- c("NA (species not present)", my.at.CEN[2:length(my.at.CEN)])
  myLabels.GF <- c("NA (species not present)", my.at.GF[2:length(my.at.GF)])
  myColorkey.CEN <- list(at=my.brks.CEN, labels=list(at=my.at.CEN, labels=myLabels.CEN))
  myColorkey.GF <- list(at=my.brks.GF, labels=list(at=my.at.GF, labels=myLabels.GF))
  myTheme <- rasterTheme(region = c(zeroCol, myPalette))
  # Plot
  levelplot(deltaCEN, par.settings = myTheme, at=my.at.CEN, colorkey=myColorkey.CEN, maxpixels = 1e8, margin = FALSE, main = "CEN diff")
  levelplot(genomicOffset, par.settings = myTheme, at=my.at.GF, colorkey=myColorkey.GF, maxpixels = 1e8, margin = FALSE, main = "GF diff")
  levelplot(glacialGenomicOffset, par.settings = myTheme, at=my.at.GF, colorkey=myColorkey.GF, maxpixels = 1e8, margin = FALSE, main = "GlacOff raw")
}

## Plot with DEM baselayer
diff_raster_stack <- stack(deltaCEN, genomicOffset, glacialGenomicOffset)
diff_dem_raster_stack <- stack()
for(i in 1:nlayers(diff_raster_stack)) {
  diff_raster.dem <- diff_raster_stack[[i]]
  diff_raster.dem[diff_raster.dem == 0.00] <- NA
  diff_raster.dem <- diff_raster.dem + 5000
  if(cast == "forecast" || i == 1 ) {
    base_map_dem_init <- GLACIER_LANDSEA_DEM_TIMESERIES[[211]]
  } else if (cast == "hindcast" && i >= 1) {
    base_map_dem_init <- GLACIER_LANDSEA_DEM_TIMESERIES[[1]]
  }
  base_map_dem_init[base_map_dem_init == -32768] <- NA # assign NA value
  diff_raster.dem <- cover(diff_raster.dem, base_map_dem_init)
  if(cast == "forecast" || i == 1 ) {
    diff_raster.dem <- overlay(diff_raster.dem, waterBodies_raster, fun=function(x, y) ifelse(is.na(y), NA, x))
  }
  proj4string(diff_raster.dem) = CRS(prj_longlat)
  # Append to stack
  diff_dem_raster_stack <- addLayer(diff_dem_raster_stack, diff_raster.dem)
}

## Define plotting parameters.
nbreaks <- 11
min_DEM_zValue <- min(minValue(GLACIER_LANDSEA_DEM_TIMESERIES)[which(minValue(GLACIER_LANDSEA_DEM_TIMESERIES) >-32768)])
max_DEM_zValue <- max(maxValue(GLACIER_LANDSEA_DEM_TIMESERIES))
breakpoints_DEM <- seq(min_DEM_zValue, max_DEM_zValue, length.out = 100)
my.at.CEN.dem <- 5000 + seq(-0.2, 0.37, length.out = nbreaks)
#my.at.GF.dem <- 5000 + seq(min(genomicOffset@data@min, glacialGenomicOffset@data@min), max(genomicOffset@data@max, glacialGenomicOffset@data@max), length.out = nbreaks)
my.at.GF.dem <- 5000 + seq(0, 0.038, length.out = nbreaks)
breakpoints_projection_CEN <- my.at.CEN.dem + ((my.at.CEN.dem[2] - my.at.CEN.dem[1])/2)
breakpoints_projection_GF <- my.at.GF.dem + ((my.at.GF.dem[2] - my.at.GF.dem[1])/2)
breakpoints_CEN <- c(breakpoints_DEM, breakpoints_projection_CEN)
breakpoints_GF <- c(breakpoints_DEM, breakpoints_projection_GF)
color.palette=colorRampPalette(c('gray95','gray50'))
custom_colors <- c(color.palette(100), rev(brewer.pal(11, "RdYlBu")))
# Plot
levelplot(diff_dem_raster_stack[[1]], col.regions = custom_colors, at=breakpoints_CEN,  maxpixels = 1e8, margin = FALSE)
levelplot(diff_dem_raster_stack[[2]], col.regions = custom_colors, at=breakpoints_GF,  maxpixels = 1e8, margin = FALSE)
levelplot(diff_dem_raster_stack[[3]], col.regions = custom_colors, at=breakpoints_GF,  maxpixels = 1e8, margin = FALSE)
