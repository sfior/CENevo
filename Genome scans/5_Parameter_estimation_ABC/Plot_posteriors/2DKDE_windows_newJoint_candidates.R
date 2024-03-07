# This script plots 2DKDE posteriors for the high-low migration rates across windows (compatible with grid and MCMC joints)

# See: https://stats.stackexchange.com/questions/191725/sample-from-distribution-given-by-histogram

# Import libraries
library(MASS)
library(RColorBrewer)
library(MASS)

# Define whether joint posteriors were calculated via grid search "grid" or MCMC sampling "MCMC"
species = "Dsyl"
joint.type = "grid"
start_window = 304
end_window = 305
#Dsyl, Chr1, CEN: 143-162
#Dsyl, Chr13, TCF1; 182-188
#Dcar, Chr8, TCF1; 1142-1149
#Dcar, Chr15, FT; 2773-2791 (ns)

# Define run parameters
density_points = 33 # number of points per parameter.
prior_migHL <- c(-3,3)
prior_migLH <- c(-3,3)

# Define neutral estimates
if (species == "Dsyl") {
  neutral_estimate <- c(3.33333, 2.57576, 2.51515, 5.82828, 5.78789, 3.57576, 3.22222, 3.0202, 2.56566, 2.81818, 2.21212) # RECON2 newSS minDP10maxDP200 3RunsCombined NEW DSYL
  neutral_estimate_2D <- c(2.57576, 2.51515) # RECON2_newSS_minDP10maxDP200_Dsyl_3RunsCombined marginalNeutral 
} else if (species == "Dcar") {
  neutral_estimate <- c(2.87879, 2.39393, 2.51515, 5.50505, 5.50508, 3.17172, 3.22222, 3.22222, 2.11112, 2.11111, 2.36364) # RECON2 newSS minDP10maxDP200 3RunsCombined NEW DCAR
  neutral_estimate_2D <- c(2.39393, 2.51515) # RECON2_newSS_minDP10maxDP200_Dcar_3RunsCombined marginalNeutral
}

# Read in ABCEstimation data
if (joint.type == "grid") {
  # For jointPosteriors
  #setwd("/Users/luqman/Desktop/GenomeScan/Estimation results/Estimation_results_simpleModel_11params_RECON2_newSS_minDP10maxDP200_Dsyl_3RunsCombined_PLS20_gW_Joint_retSims5000_gridPoints33_FINAL/jointPosteriors/")
  setwd("/Users/luqman/Desktop/GenomeScan/Estimation results/Estimation_results_simpleModel_11params_RECON2_newSSminDP10maxDP200_Dsyl_3RunsCombined_PLS20_Joint_retSims10000_gridPoints33_FINAL_Chr1/jointPosteriors/")
  if (species == "Dsyl") {
    prefix <- "ABC_estimation_simpleModel_11params_RECON2_newSSminDP10maxDP200_Dsyl_3RunsCombined_PLS20_model0_jointPosterior_2_3_Obs"
    #prefix <- "ABC_estimation_simpleModel_11params_RECON2_newSS_minDP10maxDP200_Dsyl_3RunsCombined_PLS20_gW_Joint_retSims5000_gridPoints33_FINAL_model0_jointPosterior_2_3_Obs"
  } else if (species == "Dcar") {
    prefix <- "ABC_estimation_simpleModel_11params_RECON2_newSS_minDP10maxDP200_Dcar_3RunsCombined_PLS20_model0_jointPosterior_2_3_Obs"
  }
} else if (joint.type == "MCMC") {
  # For jointPosteriorSamples (MCMC)
  setwd("/Users/luqman/Desktop/Dsyl_Chr1/jointMCMCsamples/")
  if (species == "Dsyl") {
    prefix <- "ABC_estimation_simpleModel_11params_RECON2_newSSminDP10maxDP200_Dsyl_3RunsCombined_PLS20_model0_jointPosterior_2_3_Obs"
  } else if (species == "Dcar") {
    prefix <- "ABC_estimation_simpleModel_11params_RECON2_newSS_minDP10maxDP200_Dcar_3RunsCombined_PLS20_model0_jointPosterior_2_3_Obs"
  } 
}

# Define 2D KDE plotter function. See: https://stackoverflow.com/questions/30814545/r-get-joint-probabilities-from-2d-kernel-density-estimate)
plot2DKDE <- function(joint_type, input_data, axis_min, axis_max, neutral_x, neutral_y, cred_interval, iteration) {
  # Specify whether joint (mHL and mLH) parameters were calculated via grid search or MCMC sampling
  if (joint.type == "grid") {
    # Let's calculate symmetry of posterior distribution (about the diagonal)
    # Recall that there are points above the diagonal, below the diagonal and on the diagonal. To reduce degrees of freedom, we first account for the points on the diagonal (i.e. where mLH=Mhl).
    symmetry <- sum((ABC_GLM$log_m_highcontinent_lowcontinent == ABC_GLM$log_m_lowcontinent_highcontinent)*ABC_GLM$density)
    # Then the integral of the area above the diagonal equals one minus the integral of the area below the diagonal. So we can just output one variable.
    asymmetry <- sum((ABC_GLM$log_m_highcontinent_lowcontinent < ABC_GLM$log_m_lowcontinent_highcontinent)*ABC_GLM$density) / (sum(ABC_GLM$density) - symmetry)
    # Then prepare the density plot
    dens <- list()
    # This (x-axis) is refers to log_m_highcontinent_lowcontinent
    dens[[1]] <- seq(prior_migHL[1],prior_migHL[2], length.out = density_points)
    # This (y-axis) is refers to log_m_lowcontinent_highcontinent
    dens[[2]] <- seq(prior_migLH[1],prior_migLH[2], length.out = density_points)
    # We transform the density vector into matrix.
    dens[[3]] <- matrix(ABC_GLM$density,nrow=density_points, ncol = density_points)
    names(dens) <- c("x", "y", "z")
  } else if (joint.type == "MCMC") {
    # Let's calculate symmetry of posterior distribution (about the diagonal)      
    # Recall that there are points above the diagonal, below the diagonal and on the diagonal. To reduce degrees of freedom, we first account for the points on the diagonal (i.e. where mLH=Mhl).
    symmetry <- sum(ABC_GLM$log_m_highcontinent_lowcontinent == ABC_GLM$log_m_lowcontinent_highcontinent)
    # Then the integral of the area above the diagonal equals one minus the integral of the area below the diagonal. So we can just output one variable.
    asymmetry <- sum(ABC_GLM$log_m_highcontinent_lowcontinent < ABC_GLM$log_m_lowcontinent_highcontinent) / (length(ABC_GLM$log_m_highcontinent_lowcontinent) - symmetry)
    dens <- kde2d(ABC_GLM[,2], ABC_GLM[,3], n = density_points)
  }
  # Define argument behaviour
  if(missing(axis_min) & missing(axis_max)) {
    lim <- c(min(input_data[,c(x,y)]), max(input_data[,c(x,y)]));
  } else {
    lim <- c(axis_min,axis_max)
  }
  if(missing(neutral_x) & missing(neutral_y)) {
    neutral_x <- 0
    neutral_y <- 0
  }
  if(missing(iteration)) {
    iteration <- "NA"
  }
  if(missing(cred_interval)) {
    cred_interval <- c(0.95, 0.99)
  }
  # Output intersection point of neutral estimate and 2D joint distribution
  neutral_contour_x <- which.min(abs(dens$x - neutral_estimate_2D[1]))
  neutral_contour_y <- which.min(abs(dens$y - neutral_estimate_2D[2]))
  neutral_contour_xy <- dens$z[neutral_contour_x, neutral_contour_y]
  # Calculate credible interval (here: highest density interval)
  const2d <- sum(dens$z)
  spxx2d <- sort(dens$z, decreasing = TRUE) / const2d
  crit_1 <- spxx2d[which(cumsum(spxx2d) >= cred_interval[1])[1]] * const2d
  crit_2 <- spxx2d[which(cumsum(spxx2d) >= cred_interval[2])[1]] * const2d
  # Significance
  #neutrality_inf_KDE2D <- neutral_contour_xy > crit_1
  # Additionally, we indicate at which % credible interval the neutral estimate lies (i.e. the area above the horizontal threshold defined by the neutral estimate)
  neutral_index2d <- match((neutral_contour_xy / const2d), spxx2d)
  neutral_CI_2d <- cumsum(spxx2d)[neutral_index2d]   # this variable indicates at which HDR credible interval the neutral estimate lies (in fraction)
  # Plot
  # Choose between defining a common z scale or a window-specific z-scale
  zlim <- round(0.30,2) # fix to max of the 4 windows of interest
  #zlim <- round((max(dens$z) * 1.1),2)
  # In addition to potting the contour plot, we add the neutral point, a diagonal line, and label the contour for the CIs.
  #contour(dens, xlim=lim, ylim=lim, xlab=names(input_data)[x], ylab=names(input_data)[y])
  #filled.contour(dens, xlim=lim, ylim=lim, xlab=colnames(input_data)[2], ylab=colnames(input_data)[3], cex.lab = 1, levels=seq(0, zlim, length.out = 40), 
  filled.contour(dens, xlim=lim, ylim=lim, cex.lab = 1, levels=seq(0, zlim, length.out = 40), 
                 color.palette=colorRampPalette(c('white','red3')), main=paste0("2D Joint Posterior Distribution (",joint.type,"): migHL vs migLH - Window ", i), 
                 plot.axes = {axis(1); axis(2); lines(lim, lim, lty=3, col='black'); points(neutral_x, neutral_y, pch=19, cex = 2, col="red4"); 
                   # Choose whether to add contour labels or not
                   #contour(dens, levels = c(crit_1, crit_2, neutral_contour_xy), lty = 5, lwd = 1, labcex = 0.9, labels = c((paste0(cred_interval[1]*100,"% CI: ", round(crit_1,4))),(paste0(cred_interval[2]*100,"% CI: ", round(crit_2,4))), (paste0(round(neutral_CI_2d*100),"% CI (neutral): ", round(neutral_contour_xy,4)))), col = c("deepskyblue","blue","red4"), add = TRUE)})
                   contour(dens, levels = c(crit_1, crit_2, neutral_contour_xy), lty = 5, lwd = 1, labcex = 0.9, labels = c("", "", ""), col = c("deepskyblue","blue","red4"), add = TRUE)})
  # To add line, e.g. to check symmetry. Note that when using filled.contour, x-axis scale is shifted, so we need to correct this (manually)
  # To add point, e.g. neutral parameter point estimate
  # To add relevant text information (in-box legend)
  #text(lim[1]+(lim[2]-lim[1])*0.2, lim[2]-(lim[2]-lim[1])*0.05, cex = 0.85, labels= paste0("Symmetry of x-y parameters: ", round(asymmetry, digits=4)))
  #text(lim[1]+(lim[2]-lim[1])*0.2, lim[2]-(lim[2]-lim[1])*0.1, cex = 0.85, labels = paste0("Density at neutral estimate: ", round(neutral_contour_xy,4)))
  #text(lim[1]+(lim[2]-lim[1])*0.2, lim[2]-(lim[2]-lim[1])*0.15, cex = 0.85, labels = paste0("% HDR credible interval at neutral estimate: ", round(neutral_CI_2d*100,4), "%"))
  #text(lim[1]+(lim[2]-lim[1])*0.2, lim[2]-(lim[2]-lim[1])*0.2, cex = 0.85, labels = paste0(cred_interval[1]*100,"% HDR credible interval: ", round(crit_1,4)))
  #text(lim[1]+(lim[2]-lim[1])*0.2, lim[2]-(lim[2]-lim[1])*0.25, cex = 0.85, labels = paste0("Observation diverges significantly from neutrality: ", neutral_contour_xy < crit_1))
}

# Plot 2D KDE.
par(mfrow=c(1,1),mar=c(5.1, 4.1, 4.1, 2.1))
for (i in start_window:end_window) {
  #pdf(paste0("Window",i,".pdf"), width = 10.5, height = 8)
  ABC_GLM<-read.delim(paste(prefix, i, ".txt", sep = ""))
  # Plot the 2DKDE for migHL and migLH (parameters 2 and 3)
  plot2DKDE(joint.type, ABC_GLM, -3, 3, neutral_estimate_2D[1], neutral_estimate_2D[2], c(0.95, 0.99), i)
  #dev.off()
}

