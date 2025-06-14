##### This script performs gradient forest modelling to predict compositional shifts in whole-genome adaptive genetic variation between current and future time #####

# Load libraries
library(raster)
library(gradientForest)
library(gtools)
library(FNN)
library(concaveman)
library(gtools)

# Define run parameters
working_dir <- "/Users/hl636/Documents/Hirzi/ETHPHD/PopulationGenetics_wholeGenome/Variant calling for Central Alps/"
setwd(working_dir)
# Define gradientForest run parameters
num_SNPs_perRun <- 10000
number_total_SNPs <- 390262
# Forecast model
CMIP <- 6
# Define plotting parameters
plot_figures <- "PARTIAL" # ALL, PARTIAL, NONE
## Import population coordinate data
pop_geodata <- read.csv(paste0(working_dir,"popInfo_w14indsCorrectHe_NEW_RevisedR2ENV.csv"))
pop_geodata_GF <- pop_geodata
# Define whether object is of class gradientForest or combinedGradientForest, i.e. whether you want to run gradientForest (FALSE) or whether you want to use an imported a combinedGradientForest object (TRUE)
combined_gradientForest <- TRUE
prj_longlat <- "+init=epsg:4326"
e_alps <- extent(5, 12, 43.5, 48)
#ELEV.alps <- crop(ELEV.alps, e_alps)
GLACIER_LANDSEA_DEM_TIMESERIES <-stack("/Users/hl636/Documents/Hirzi/ETHPHD/SDM/Dsylvestris/Data/Chelsa_TraCE_Glacier_LandSea_DEM_TimeSeries.grd")
GLACIER_LANDSEA_DEM_TIMESERIES <- crop(GLACIER_LANDSEA_DEM_TIMESERIES, e_alps)
ENV_RASTERS_forWater <-stack("./ENV_RASTERS_TODAY_EUROPE_CHELSA_LARGE.grd") 
waterBodies_raster <- ENV_RASTERS_forWater$PH_5cm
waterBodies_raster <- crop(waterBodies_raster, e_alps)
ELEV.alps <- GLACIER_LANDSEA_DEM_TIMESERIES[[211]]
rm(GLACIER_LANDSEA_DEM_TIMESERIES)
ELEV.alps[ELEV.alps == -32768] <- NA # assign NA value
crs(ELEV.alps) <- crs(prj_longlat)
ELEV.alps <- overlay(ELEV.alps, waterBodies_raster, fun=function(x, y) ifelse(is.na(y), NA, x))

# If running gradientForest:
if (combined_gradientForest == FALSE) { 
  ## Import allele frequency data
  allele_data <- read.table(paste0(working_dir,"AllPops_innerMerged.minQ30minDP1204minAC14miss10_moreStringent.prim.SNP.exons_minDepth7.AF.AFcolumnsOnly"))
  # Transpose
  allele_data <- as.data.frame(t(allele_data))
  # Convert first row to header
  names(allele_data) <- as.character(unlist(allele_data[1,]))
  allele_data <- allele_data[-1,]
  colnames(allele_data)[1] <- "pop"
  # Rename pops (remove "AF_" prefix)
  allele_data$pop <- gsub('AF_', '', allele_data$pop)
  
  ## Merge population coordinate and SNP allele frequency data
  allele_data <- merge(pop_geodata, allele_data, by.x="pop", by.y="pop", all.x=TRUE)
  # Convert first column to rownames
  allele_data_temp <- allele_data[,-1]
  rownames(allele_data_temp) <- allele_data[,1]
  allele_data <- allele_data_temp
  rm(allele_data_temp)
  # Convert all classes to numeric
  allele_data[] <- lapply(allele_data, function(x) {
    if(is.factor(x)) as.numeric(as.character(x)) else x
  })
  #sapply(allele_data, class)
}

## Import current environmental data
CLIM.alps <-stack(paste0(working_dir,"ENV_RASTERS_TODAY_EUROPE_CHELSA_LARGE.grd"))
# Choose predictor variables. 
# Here we take uncorrelated variables (< 0.7) assessed at a larger extent of occurrences (SDM occurrences) and evaluated using pairs.hist.cor(Occurrences.Climate[, 4:28], cor.method = "pearson"), and round(cor(Occurrences.Climate[, 4:28]),2); varCor <- cor(Occurrences.Climate[, 4:28]); allDistNew <- abs(as.dist(cor(Occurrences.Climate[, 4:28]))); allClusNew <- hclust(1 - allDistNew); plot(allClusNew, hang=-1)
# To choose among the correlated variables, we refer to the importance of each variable evaluated from full (all 27 predictor variables) GDM and GF models (preliminary analyses). In case there is disagreement, we choose the optimum (most informative) compromise. For more details, see Predictor variable selection.xlsx
# We also assess the variance inflation factors (VIF) of the variables (conditional on VIF < 10).
pred.var.selection <- c("TEMP_DIURNAL_RANGE", "TEMP_SEASONALITY", "TEMP_MIN_COLDEST_MONTH", "TEMP_MEAN_WETTEST_QUARTER", "TEMP_MEAN_DRIEST_QUARTER", "PREC_SEASONALITY", "PREC_WARMEST_QUARTER", "PREC_COLDEST_QUARTER", "PH_5cm", "SLOPE", "LONGITUDE", "LATITUDE")
CLIM.alps <- subset(CLIM.alps, pred.var.selection)
# And crop to desired geographic extent
e <- extent(4, 15, 43, 48)
CLIM.alps <- crop(CLIM.alps, e)
CLIM.alps.gf <- CLIM.alps
# Similarly, import future environmental data
if (CMIP == 5) {
  CLIM.alps.future <-stack("/Users/hl636/Documents/Hirzi/ETHPHD/SDM/Dsylvestris/Data/ENV_RASTERS_CMIP5_EUROPE_ENSEMBLE_ALL_LARGE.grd")
} else if (CMIP == 6) {
  CLIM.alps.future <-stack("/Users/hl636/Documents/Hirzi/ETHPHD/SDM/Dsylvestris/Data/ENV_RASTERS_CMIP6_EUROPE_ENSEMBLE_ALL_LARGE.grd")
}
CLIM.alps.future <- subset(CLIM.alps.future, pred.var.selection)
CLIM.alps.future <- crop(CLIM.alps.future, e)

if (combined_gradientForest == FALSE) {
  
  ## Combine environmental data and allele frequency data
  CLIM.extract <- data.frame(extract(CLIM.alps, allele_data[,c("x","y")]))
  CLIM_allele <- data.frame(CLIM.extract, allele_data) # let's use lat, lon as a proxy for PCA1, PCA2 (as correlation between these are high)
  preds <- c(names(CLIM.extract))
  loci <- names(allele_data[,-c(1,2)])

# To speed up analysis, we run gradientForest in parallel subsets (e.g. on the cluster)
  idx_start <- ((run_number-1)*num_SNPs_perRun) + 1
  idx_end <- run_number*num_SNPs_perRun
  if (run_number < ceiling(number_total_SNPs/num_SNPs_perRun)) {
    loci <- loci[idx_start:idx_end]
  } else if (run_number == ceiling(number_total_SNPs/num_SNPs_perRun)) {
    loci <- loci[idx_start:number_total_SNPs]
  }

  ## Run gradientForest
  gfVars <- gradientForest(CLIM_allele,
                           predictor.vars=preds,
                           response.vars=loci,
                           ntree=500,
                           maxLevel=floor(log2(0.368*nrow(CLIM_allele)/2)),
                           trace=T,
                           corr.threshold=0.50) # or 0.7
  # Often a warning will prompt. Type 'warnings()' to see them. Warning refers to whether you want to use regression on columns with fewer than n unique values.
  
  # Save output
  #time_stamp <- gsub(":","_",gsub(" ", "_", Sys.time()))"
  #saveRDS(gfVars, paste0("gfModel_selectedPredVars.split", run_number,".rds"))
  #gfVars <- readRDS(paste0("gfModel_allPredVars.",time_stamp,".rds"))
}


# ## Combine parallel (subset) runs of gradientForest into one gradientForest model (object).
# # This function is very resource heavy (>100GB memory for this job) so run on Euler!
#  
# # List gradientForest parallel runs
# #files <- list.files(path = working_dir, pattern = "^gfModel_allPredVars.split.*\\.10000SNPsets.rds$")
# files <- list.files(path = working_dir, pattern = "^gfModel_selectedPredVars.split.*\\.rds$")
# files <- mixedsort(files)
# 
# # Read results into R
# gfModel_subset_list <- lapply(files, readRDS)
# 
# # Combine parallel (subset) runs into an aggregate gradientForest model (object)
# gfModel_combined <- combinedGradientForest(gfModel_subset_list[[1]],gfModel_subset_list[[2]],gfModel_subset_list[[3]],gfModel_subset_list[[4]],gfModel_subset_list[[5]],gfModel_subset_list[[6]],gfModel_subset_list[[7]],gfModel_subset_list[[8]],gfModel_subset_list[[9]],gfModel_subset_list[[10]],gfModel_subset_list[[11]],gfModel_subset_list[[12]],gfModel_subset_list[[13]],gfModel_subset_list[[14]],gfModel_subset_list[[15]],gfModel_subset_list[[16]],gfModel_subset_list[[17]],gfModel_subset_list[[18]],gfModel_subset_list[[19]],gfModel_subset_list[[20]],gfModel_subset_list[[21]],gfModel_subset_list[[22]],gfModel_subset_list[[23]],gfModel_subset_list[[24]],gfModel_subset_list[[25]],gfModel_subset_list[[26]],gfModel_subset_list[[27]],gfModel_subset_list[[28]],gfModel_subset_list[[29]],gfModel_subset_list[[30]],gfModel_subset_list[[31]],gfModel_subset_list[[32]],gfModel_subset_list[[33]],gfModel_subset_list[[34]],gfModel_subset_list[[35]],gfModel_subset_list[[36]],gfModel_subset_list[[37]],gfModel_subset_list[[38]],gfModel_subset_list[[39]],gfModel_subset_list[[40]])
# 
# # Save object
# #saveRDS(gfModel_combined, "gfModel_allPredVars_combined_allSplits.rds")
# saveRDS(gfModel_combined, "gfModel_selectedPredVars_combined_allSplits.rds")

#if (combined_gradientForest == TRUE) { 
# Choose default options (method=2, standardize="before"). Choose between cor=0.5,0.7 - results appear very similar.
gfVars <- readRDS(paste0(working_dir,"gfModel_cor0.5_selectedPredVars_woSolarRad_combined_method2standardizeBefore_allSplits.rds"))

# Print results
# Predictor importance (recall that the accuracy important is how well each variable is at predicting the response)
vars <- names(importance(gfVars))
if (plot_figures == "ALL") {
  importance(gfVars)
  plot(gfVars, plot.type="O")
}

# Make df of variable importance
var_importance <- gfVars$imp.rsq
var_importance <- var_importance[,-1]
var_importance.df <- as.data.frame(rowSums(var_importance))
var_importance.df[,2] <- 0
#var_importance.sorted.df <- var_importance.df[rev(order(var_importance.df[1])),]
var_importance.sorted.df <- var_importance.df[rev(order(var_importance.df[,1])),]
predVar_colSum <- colSums(var_importance.sorted.df)[1]
var_importance.sorted.df[,2] <- var_importance.sorted.df[,1]/predVar_colSum
colnames(var_importance.sorted.df) <- c("importance", "normalized_importance")

if (plot_figures == "ALL") {
  # Measure of model fit (R2) by SNP (gradientForest) or SNP set (combinedGradientForest)
  if (combined_gradientForest == FALSE) {
    plot(gfVars, plot.type="P", show.names=T, horizontal=F)
  } else if (combined_gradientForest == TRUE) {
    #plot(gfVars, plot.type="Performance") # default plot
    # Extract total-R2 and ENV(only)-R2 from gfvars. The former can be found directly at gfVars$rsq, however the latter needs to be calculated.
    r2_allpreds <- as.data.frame(t(gfVars$imp.rsq))[-c(1),]
    r2_allpreds[] <- lapply(r2_allpreds, function(x) {
      if(is.factor(x)) as.numeric(as.character(x)) else x
    })
    #sapply(r2_allpreds, class)
    # Calculate total-R2 and ENV(only)-R2
    r2_allpreds[,ncol(r2_allpreds)+1] <- rowSums(r2_allpreds[,1:10])
    colnames(r2_allpreds)[ncol(r2_allpreds)] <- "r2_ENV"
    r2_allpreds[,ncol(r2_allpreds)+1] <- rowSums(r2_allpreds[,1:12])
    colnames(r2_allpreds)[ncol(r2_allpreds)] <- "r2_total"
    # Format data for boxplot
    #r2_factors <- split(gfVars$rsq,factor(rep(names(gfVars$nspec),gfVars$nspec)))
    #r2_factors <- split(r2_allpreds$r2_total,factor(rep(names(gfVars$nspec),gfVars$nspec))) # for total R2, identical to above
    r2_factors <- split(r2_allpreds$r2_ENV,factor(rep(names(gfVars$nspec),gfVars$nspec))) # for env-only R2 (i.e. excluding latlon)
    r2_factors.sorted <- r2_factors[rev(mixedsort(names(r2_factors)))] # for vertical
    boxplot(r2_factors.sorted,plot=TRUE, col=rgb(0.5,0.65,0.85,0.85), pars=list(outpch = 21, outcol="grey10", outbg = "brown1", outcex=1.65), las=2, xlab="R2", ylab="gradientForest runs", horizontal = TRUE, ylim = c(0, 0.9)) # nicer plot
  }
  
  # Plot split densities of predictor variables
  if (combined_gradientForest == FALSE) {
    plot(gfVars, plot.type="S")
  } else if (combined_gradientForest == TRUE) {
    plot(gfVars, plot.type="Predictor.Density")
  }

  # Cumulative importance plot
  if (combined_gradientForest == FALSE) {
    # Species-level  
    plot(gfVars, plot.type="C", show.overall=F, leg.nspecies=5, cex.lab=1, common.scale=T, par.args=list(mgp=c(1.5,0.5,0), mar=c(2.5,1,0.1,0.5), omi=c(0,0.3,0,0)))
    # Community-level
    plot(gfVars, plot.type="C", show.species=F, leg.nspecies=5, cex.lab=1, common.scale=T, par.args=list(mgp=c(1.5,0.5,0), mar=c(2.5,1,0.1,0.5), omi=c(0,0.3,0,0)))
  } else if (combined_gradientForest == TRUE) {
    # Cumulative importance plot - species-level  
    plot(gfVars, plot.type="C", cex.lab=1, common.scale=T, par.args=list(mgp=c(1.5,0.5,0), mar=c(2.5,1,0.1,0.5), omi=c(0,0.3,0,0)))
  }
}

## Project model in space
CLIM.alps.df <- as.data.frame(CLIM.alps, xy=TRUE, na.rm=TRUE)
CLIM.alps.df.GF <- data.frame(CLIM.alps.df[,c("x","y")], predict(gfVars,CLIM.alps.df[,vars]))

# The multi-dimensional biological space can most effectively be represented by taking the principle components of the transformed grid and presenting the first two dimensions in a biplot.
# To remove neutral structure (here represented in proxy via longitude and latitude) from the visualisation, we calculate and plot the PCs excluding contributions from latitude and longitude.
vars.woLatLon <- vars [! vars %in% c("LATITUDE", "LONGITUDE")]
vars.woLatLon.gf <- vars.woLatLon
# Consider if you want to scale the PCA, or if you'd like to retain/preserve differences in the magnitude of genetic importance among the environmental variables (Fitzpatrick & Keller 2015, pg.6). Let's choose the latter (without scaling).
PCs <- prcomp(CLIM.alps.df.GF[, vars.woLatLon], scale. = FALSE, center = TRUE)
PCs.full.gf <- PCs # Output for ensemble model generation
PCs.gf <- PCs$x # Output for ensemble model generation

# Let's convert this to a raster
PCs.df <- as.data.frame(PCs$x)
PCs.df <- cbind(CLIM.alps.df.GF[,c("x","y")], PCs.df)
PCs.rast <- rasterFromXYZ(PCs.df)
crs(PCs.rast) <- crs(CLIM.alps)

## Forecast model
# Project model in future space
CLIM.alps.future.df <- as.data.frame(CLIM.alps.future, xy=TRUE, na.rm=TRUE)
CLIM.alps.future.df.GF <- data.frame(CLIM.alps.future.df[,c("x","y")], predict(gfVars,CLIM.alps.future.df[,vars]))
# Then calculate the PCA of the future biological space. We do so in the PCA space of the current climate (i.e. that this future PCA will be transformed, scaled and centered based on the current climate PCA), for apples-to-apples comparison
CLIM.alps.future.df.GF.tranformed <- scale(CLIM.alps.future.df.GF[, rownames(PCs$rotation)], scale = FALSE, center = PCs$center)
# Consider if you want to scale the PCA, or if you'd like to retain/preserve differences in the magnitude of genetic importance among the environmental variables (Fitzpatrick & Keller 2015, pg.6). Let's choose the latter (without scaling).
#CLIM.alps.future.df.GF.tranformed <- scale(CLIM.alps.future.df.GF.tranformed, scale = PCs$scale, center = FALSE)
PCs.future <-  as.matrix(CLIM.alps.future.df.GF.tranformed) %*% as.matrix(PCs$rotation)
PCs.future.gf <- PCs.future # Output for ensemble model generation

# Let's convert this to a raster
PCs.future.df <- as.data.frame(PCs.future)
PCs.future.df <- cbind(CLIM.alps.future.df.GF[,c("x","y")], PCs.future.df)
PCs.future.rast <- rasterFromXYZ(PCs.future.df)
crs(PCs.future.rast) <- crs(CLIM.alps.future)

# Define colour palette - scale rasters by common scale
# Current
r.scaled <- (PCs.rast[[1]]-min(PCs.rast[[1]]@data@min, PCs.future.rast[[1]]@data@min)) / (max(PCs.rast[[1]]@data@max, PCs.future.rast[[1]]@data@max)-min(PCs.rast[[1]]@data@min, PCs.future.rast[[1]]@data@min))*255
g.scaled <- (PCs.rast[[2]]-min(PCs.rast[[2]]@data@min, PCs.future.rast[[2]]@data@min)) / (max(PCs.rast[[2]]@data@max, PCs.future.rast[[2]]@data@max)-min(PCs.rast[[2]]@data@min, PCs.future.rast[[2]]@data@min))*255
b.scaled <- (PCs.rast[[3]]-min(PCs.rast[[3]]@data@min, PCs.future.rast[[3]]@data@min)) / (max(PCs.rast[[3]]@data@max, PCs.future.rast[[3]]@data@max)-min(PCs.rast[[3]]@data@min, PCs.future.rast[[3]]@data@min))*255
pcaRast.rgb <- stack(r.scaled, g.scaled, b.scaled)
# future
r.future.scaled <- (PCs.future.rast[[1]]-min(PCs.rast[[1]]@data@min, PCs.future.rast[[1]]@data@min)) / (max(PCs.rast[[1]]@data@max, PCs.future.rast[[1]]@data@max)-min(PCs.rast[[1]]@data@min, PCs.future.rast[[1]]@data@min))*255
g.future.scaled <- (PCs.future.rast[[2]]-min(PCs.rast[[2]]@data@min, PCs.future.rast[[2]]@data@min)) / (max(PCs.rast[[2]]@data@max, PCs.future.rast[[2]]@data@max)-min(PCs.rast[[2]]@data@min, PCs.future.rast[[2]]@data@min))*255
b.future.scaled <- (PCs.future.rast[[3]]-min(PCs.rast[[3]]@data@min, PCs.future.rast[[3]]@data@min)) / (max(PCs.rast[[3]]@data@max, PCs.future.rast[[3]]@data@max)-min(PCs.rast[[3]]@data@min, PCs.future.rast[[3]]@data@min))*255
pcafutureRast.rgb <- stack(r.future.scaled, g.future.scaled, b.future.scaled)

## Define colour palette
color.palette=colorRampPalette(c('gray95','gray50'))

## Map in geographic space
crs(pcaRast.rgb) <- crs(CLIM.alps)
crs(pcafutureRast.rgb) <- crs(CLIM.alps.future)
# Crop
pcaRast.rgb <-  crop(pcaRast.rgb, e_alps)
pcafutureRast.rgb <-  crop(pcafutureRast.rgb, e_alps)
# The following map plots predicted PC scores in geographic coordinates and represents continuous changes in inferred compositional patterns associated with the predictors.
if (plot_figures == "ALL" || plot_figures == "PARTIAL" ) {
  par(mfrow=c(1,1), mar = c(4, 4, 3, 3))
  plotRGB(pcaRast.rgb, r=1, g=2, b=3)
  plotRGB(pcafutureRast.rgb, r=1, g=2, b=3)
}

# ## Biplot of biological space
# # Define colour palette
# r.biplot <- PCs$x[, 1]
# g.biplot <- PCs$x[, 2]
# b.biplot <- PCs$x[, 3]
# # Scale between 0-255. Use same standardisation scale as above.
# r.biplot <- (r.biplot - min(PCs.rast[[1]]@data@min, PCs.future.rast[[1]]@data@min)) / (max(PCs.rast[[1]]@data@max, PCs.future.rast[[1]]@data@max)-min(PCs.rast[[1]]@data@min, PCs.future.rast[[1]]@data@min))*255
# g.biplot <- (g.biplot - min(PCs.rast[[2]]@data@min, PCs.future.rast[[2]]@data@min)) / (max(PCs.rast[[2]]@data@max, PCs.future.rast[[2]]@data@max)-min(PCs.rast[[2]]@data@min, PCs.future.rast[[2]]@data@min))*255
# b.biplot <- (b.biplot - min(PCs.rast[[3]]@data@min, PCs.future.rast[[3]]@data@min)) / (max(PCs.rast[[3]]@data@max, PCs.future.rast[[3]]@data@max)-min(PCs.rast[[3]]@data@min, PCs.future.rast[[3]]@data@min))*255
# 
# # The environmental variables may be shown as vectors, e.g. here limited to the most important predictors — in this example, variables to show as vectors are selected.
# nvs <- dim(PCs$rotation)[1]
# lv <- length(vars.woLatLon)
# vind <- rownames(PCs$rotation) %in% vars.woLatLon
# scal <- 40
# xrng <- range(PCs$x[, 1], PCs$rotation[, 1]/scal) * 1.1
# yrng <- range(PCs$x[, 2], PCs$rotation[, 2]/scal) * 1.1
# # And plot. Different coordinate positions in the biplot represent differing compositions, as associated with the predictors.
# par(mar = c(2, 2, 2, 2), mfrow = c(1, 1), bg=NA)
# plot((PCs$x[, 1:2]), xlim = xrng, ylim = yrng, pch = ".", cex = 4, col = rgb(r.biplot, g.biplot, b.biplot, max = 255), asp = 1, xaxt="n", yaxt="n")
# points(PCs$rotation[!vind, 1:2]/scal, pch = "+")
# arrows(rep(0, lv), rep(0, lv), PCs$rotation[vars.woLatLon, 1]/scal, PCs$rotation[vars.woLatLon, 2]/scal, length = 0.1, lwd = 2.5)
# jit <- 0.0015
# text(PCs$rotation[vars.woLatLon, 1]/scal + jit * sign(PCs$rotation[vars.woLatLon, 1]), PCs$rotation[vars.woLatLon, 2]/scal + jit * sign(PCs$rotation[vars.woLatLon, 2]), labels = vars.woLatLon)
# # Draw axes through origin
# axis(1, pos=0)
# axis(2, pos=0)
# # Add the location of populations in biological space (circles)
# # First transform the population environmental predictors, which are available from gf$X.
# Trns_site <- predict(gfVars)
# PCsites <- predict(PCs, Trns_site[, vars.woLatLon])
# PCsites <- unique(PCsites) # since we combined n gradientForest objects via combinedGradientForest, we have n*pops; i.e. duplicate rows.
# points(PCsites[, 1:2])
# # Add the weight mean location of species (crosses).
# if (combined_gradientForest == FALSE) {
#   SpsWtd <- sweep(gfVars$Y, 2, apply(gfVars$Y, 2, min), "-")
#   SpsWtdPCs <- (t(SpsWtd) %*% (PCsites[, 1:2]))/colSums(SpsWtd)
#   points(SpsWtdPCs, col = "red", pch = "+")
#   #If required, the frequency of any given SNP may be plotted on the biplot. For example the first SNP from gf$Y:
#   #sp <- colnames(SpsWtd)[1]
#   #points(PCsites[, 1:2], col = "blue", cex = SpsWtd[, sp]/2)
# }
