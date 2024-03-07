# This script plots the combined and aggregate (superimposed) marginal posterior densities for the demographic parameters for model: simpleModel_11params

# Import libraries
library(MASS)
library(RColorBrewer)

# Define input variables
species <- "Dcar" # Dsyl, Dcar or Pseudo
rescale <- TRUE # For final plotting, set rescale to true to rescale Mic to mean island Ne (the recipient)
nparams <- 11 # number of parameters
# Define plot name labels
plot_labels <- c("log10(mIC)", "log10(mHL)", "log10(mLH)", "log10(nH)", "log10(nL)", "log10(nL1)", "log10(nL2)", "log10(nL3)", "log10(nH1)", "log10(nH2)", "log10(nH3)")

if (species == "Dsyl") {
  # Number of loci; here 0-indexed so no of loci - 1
  num_loci <- 1007
  # Directory and prefix
  setwd("/Users/luqman/Desktop/Estimation_results_simpleModel_11params_RECON2_newSS_minDP10maxDP200_Dsyl_3RunsCombined_PLS20_gW_Joint_retSims5000_gridPoints33_FINAL/marginalPosteriorDensities/")
  prefix <- "ABC_estimation_simpleModel_11params_RECON2_newSS_minDP10maxDP200_Dsyl_3RunsCombined_PLS20_gW_Joint_retSims5000_gridPoints33_FINAL_model0_MarginalPosteriorDensities_Obs"
} else if (species == "Dcar") {
  num_loci <- 1049
  setwd("/Users/luqman/Desktop/Estimation_results_simpleModel_11params_RECON2_newSS_minDP10maxDP200_Dcar_3RunsCombined_PLS20_gW_Joint_retSims5000_gridPoints33_FINAL/marginalPosteriorDensities/")
  prefix <- "ABC_estimation_simpleModel_11params_RECON2_newSS_minDP10maxDP200_Dcar_3RunsCombined_PLS20_gW_Joint_retSims5000_gridPoints33_FINAL_model0_MarginalPosteriorDensities_Obs"
} else if (species == "Pseudo") {
  num_loci <- 999 
  setwd("/Users/luqman/Desktop/ABC GS Revisions/Estimation_results_simpleModel_11params_RECON2_newSS2_minDP10maxDP200_Dsyl_PLS20_revisedJoint_gridPoints33_retSims3000_PseudoObs/marginalPosteriorDensities/") #Dsyl newSS2 RECON2 pseudoObs
  prefix <- "ABC_estimation_simpleModel_11params_RECON2_newSS2_minDP10maxDP200_Dsyl_PLS20_model0_MarginalPosteriorDensities_Obs"   #Dsyl newSS2 RECON2 pseudoObs
}

####### Plot combined posteriors and all observations #######

# Find product of probability densities (take log)
prod <- matrix(0, ncol=nparams, nrow=100)
for (i in 0:num_loci) {
  ABC_GLM<-read.delim(paste(prefix, i, ".txt", sep = "")) 
  for(p in 1:nparams){
    prod[,p] <- prod[,p] + log(ABC_GLM[,2*p+1]);
  }
}

# Normalise product of probability densities
for(p in 1:nparams){
  # Normalise and plot product of densities
  prod[,p] <- prod[,p] - max(prod[,p]);
  prod[,p] <- exp(prod[,p])
  #prod[,p] <- prod[,p] / (num_loci)
}

# Get values at peaks of product of probability densities (values)
combined_estimate_num <- vector()
for(p in 1:nparams){
  comb_estimate_param_num <- (ABC_GLM[,2*p][match(max(exp(prod[,p])),exp(prod[,p]))])
  combined_estimate_num[[p]] <- comb_estimate_param_num
}

if (rescale == TRUE) {
  # For plotting, rescale Mic to mean island Ne (the recipient)
  mean_mIC_N <- mean(combined_estimate_num[6:11])
  mIC_scaling_factor <- mean_mIC_N - 4 # from RECON2 def: m_island_continent = pow10(log_m_island_continent); to rescaled (correct) def: m_island_continent = pow10(log_m_island_continent - log_N_island + 4)
  ABC_GLM[,2] <- ABC_GLM[,2] - mIC_scaling_factor
} else if (rescale == FALSE) {
  mIC_scaling_factor <- 0
}

# Get names and values at peaks of product of probability densities (names and values)
combined_estimate <- list()
for(p in 1:nparams){
  #comb_estimate_param <- (paste(colnames(ABC_GLM)[2*p], round(ABC_GLM[,2*p][match(max(exp(prod[,p])),exp(prod[,p]))],3), sep = ": "))
  comb_estimate_param <- (paste(plot_labels[p], round(ABC_GLM[,2*p][match(max(exp(prod[,p])),exp(prod[,p]))],5), sep = ": "))
  combined_estimate[[p]] <- comb_estimate_param
}

# Plot all together
par(mfrow=c(2,3))
for(p in 1:5){
  # Normalise and plot product of densities
  #plot(ABC_GLM[,2*p], prod[,p], type='l', lty=2, col='darkred', main=names(ABC_GLM[2*p]), ylim = c(0,(max(ABC_GLM[[2*p + 1]])*2)), lwd = 1, xlab = combined_estimate[[p]]);
  plot(ABC_GLM[,2*p], prod[,p], type='l', lty=2, col='darkred', ylim = c(0,1), lwd = 1, xlab = combined_estimate[[p]], ylab = NA, cex.lab = 2, cex.axis = 2);
  # Plot all windows
  lines(ABC_GLM[[2*p]],ABC_GLM[[2*p + 1]],type='l',col='red', lwd = 0.05)
  for (i in 0:num_loci) {
    ABC_GLM<-read.delim(paste(prefix, i, ".txt", sep = "")) 
    ABC_GLM[,2] <- ABC_GLM[,2] - mIC_scaling_factor # rescale mIC
    lines(ABC_GLM[[2*p]],ABC_GLM[[2*p + 1]],type='l',col='red', lwd = 0.05)
  }
}

par(mfrow=c(2,3))
for(p in 6:nparams){
  # Normalise and plot product of densities
  #plot(ABC_GLM[,2*p], prod[,p], type='l', lty=2, col='darkred', main=names(ABC_GLM[2*p]), ylim = c(0,(max(ABC_GLM[[2*p + 1]])*2)), lwd = 1, xlab = combined_estimate[[p]]);
  plot(ABC_GLM[,2*p], prod[,p], type='l', lty=2, col='darkred', ylim = c(0,0.6), lwd = 1, xlab = combined_estimate[[p]], ylab = NA, cex.lab = 2, cex.axis = 2);
  # Plot all windows
  lines(ABC_GLM[[2*p]],ABC_GLM[[2*p + 1]],type='l',col='red', lwd = 0.05)
  for (i in 0:num_loci) {
    ABC_GLM<-read.delim(paste(prefix, i, ".txt", sep = "")) 
    lines(ABC_GLM[[2*p]],ABC_GLM[[2*p + 1]],type='l',col='red', lwd = 0.05)
  }
}

print(combined_estimate)
