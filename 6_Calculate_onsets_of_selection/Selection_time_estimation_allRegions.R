
#### This script plots the selection time estimates (migTime models) for multiple candidate regions.
## Note 1: It takes as input the output of ABCestimator (normal and indrepMCMC)
## Note 2:  It outputs plots of mHL, mLH and t_mig, as well as writes a summary (results) csv file which can be used for plotting (via Plot_combined_time_posterior.R) 

### Define global input variables
# These define the candidate regions analysed
#candidate_regions <- c("DSYL_CEN", "DSYL_TCF1", "DCAR_TCF1", "DCAR_FT")
candidate_regions <- c("CEN_DSYL", "TCF1_DSYL", "TCF1_DCAR", "FT_DCAR")
# These are the number of (5000bp) windows that define a candidate region
candidate_regions_numWindows <- c(5,2,3,1)

### Run (loop) across candidate regions
for (k in 1:4) {
    
  ### Define local variables and parameters of interest
  # Set working directory (for marginal posterior probabilities files)
  #setwd(paste0("/Users/luqman/Desktop/ABCtoolbox/Estimation_results_simpleModel_11params_RECON2_migTime_FINAL_RUN2_FIXED_",candidate_regions[k])) # migTime FINAL RUN2 FIXED
  #setwd(paste0("/Users/luqman/Desktop/ABCtoolbox/Estimation_results_simpleModel_11params_RECON2_migTime_FINAL_RUN2_FREE_",candidate_regions[k])) # migTime FINAL RUN2 FREE
  setwd(paste0("/Users/luqman/Desktop/Migration time demography/Estimation_results_simpleModel_11params_RECON2_migTime_FINAL_RUN2_newSS_FREE_",candidate_regions[k])) # migTime FINAL RUN2 newSS FREE
  
  #prefix <- "ABC_estimation_simpleModel_11params_RECON2_migTime_FINAL_RUN2_FIXED_PLS10_model0_MarginalPosteriorDensities_Obs" # migTime FINAL RUN2 FIXED
  #prefix <- "ABC_estimation_simpleModel_11params_RECON2_migTime_FINAL_RUN2_FREE_PLS20_model0_MarginalPosteriorDensities_Obs" # migTime FINAL RUN2 FREE
  prefix <- "ABC_estimation_simpleModel_11params_RECON2_migTime_FINAL_RUN2_newSS_FREE_PLS20_model0_MarginalPosteriorDensities_Obs" # migTime FINAL RUN2 newSS FREE
  
  # MCMC file
  #mcmc_input <- paste0("/Users/luqman/Desktop/ABCtoolbox/Estimation_results_simpleModel_11params_RECON2_migTime_FINAL_RUN2_FIXED_",candidate_regions[k],"_indrepMCMC/ABC_estimation_simpleModel_11params_RECON2_migTime_FINAL_RUN2_FIXED_PLS10_model0_MCMC.txt") # migTime FINAL RUN2 FIXED
  #mcmc_input <- paste0("/Users/luqman/Desktop/ABCtoolbox/Estimation_results_simpleModel_11params_RECON2_migTime_FINAL_RUN2_FREE_",candidate_regions[k],"_indrepMCMC/ABC_estimation_simpleModel_11params_RECON2_migTime_FINAL_RUN2_FREE_PLS20_model0_MCMC.txt") # migTime FINAL RUN2 FREE
  mcmc_input <- paste0("/Users/luqman/Desktop/Migration time demography/Estimation_results_simpleModel_11params_RECON2_migTime_FINAL_RUN2_newSS_FREE_",candidate_regions[k],"_indrepMCMC/ABC_estimation_simpleModel_11params_RECON2_migTime_FINAL_RUN2_newSS_FREE_PLS20_model0_MCMC.txt") # migTime FINAL RUN2 newSS FREE
  
  # Model parameters
  num_loci <- candidate_regions_numWindows[k]
  nparams <- 3
  
  ### Calculate combined posteriors (product of densities method)
  prod_1D <- matrix(0, ncol=nparams, nrow=100)
  for (i in 0:(num_loci-1)) {
    ABC_GLM<-read.delim(paste(prefix, i, ".txt", sep = "")) 
    #ABC_rej<-read.delim(paste(prefix2, i, ".txt", sep = ""))  
    # For product of densities.
    for(p in 1:nparams){
      prod_1D[,p] <- prod_1D[,p] + log(ABC_GLM[,2*p+1])
    }
  }
  
  ### Import combined posteriors (MCMC independent replicates method)
  mcmc <- read.table(mcmc_input, header=TRUE)
  
  ### Plot combined posteriors (indrep and product of densities methods) & all windows in one plot
  par(mfrow=c(1,3))
  for(p in 1:nparams){
    # Normalise and plot product of densities
    prod_1D[,p] <- prod_1D[,p] - max(prod_1D[,p]);
    prod_1D[,p] <- exp(prod_1D[,p])
    #plot(ABC_GLM[,2*p], prod_1D[,p], type='l', col='green', main=names(ABC_GLM[2*p]), lwd = 2, lty = 1, xlab = NA, ylab = NA, ylim = c(0,0.5));
    plot(ABC_GLM[,2*p], prod_1D[,p], type='l', col='green', main=paste0(candidate_regions[k], "; ", names(ABC_GLM[2*p])), lwd = 2, lty = 1, xlab = NA, ylab = NA);
    # Plot indrep MCMC
    dens.indRepMCMC <- density(mcmc[,p]);
    dens.indRepMCMC$y <- dens.indRepMCMC$y / max(dens.indRepMCMC$y);
    lines(dens.indRepMCMC$x, dens.indRepMCMC$y, col='blue', lwd = 2, lty = 1)
    # Plot all windows
    lines(ABC_GLM[[2*p]],ABC_GLM[[2*p + 1]],type='l',col='red', lwd = 0.75, main=names(ABC_GLM[2*p]), ylim = c(0,(max(ABC_GLM[[2*p + 1]])*2.5)))
    for (i in 0:(num_loci-1)) {
      ABC_GLM<-read.delim(paste(prefix, i, ".txt", sep = ""))
      lines(ABC_GLM[[2*p]],ABC_GLM[[2*p + 1]],type='l',col='red', lwd = 0.75)
    }
    legend("topleft", legend=c("Window posteriors", "Combined posterior (indrepMCMC)", "Combined posterior (prod.densities)"),
           col=c("red", "blue","green"), lty=c(1,1,1), lwd = c(0.75, 1, 1), cex=0.8)
    # Print peak
    print(paste0(colnames(ABC_GLM)[2*p], "; prodOfDensities: ", ABC_GLM[,2*p][match(max(exp(prod_1D[,p])),exp(prod_1D[,p]))]))
    print(paste0(colnames(ABC_GLM)[2*p], "; indrepMCMC: ", dens.indRepMCMC$x[match(max(dens.indRepMCMC$y),dens.indRepMCMC$y)]))
  }
  
  ### Append indrep MCMC results to a dataframe
  if (k == 1) {
    selection.df <- as.data.frame(dens.indRepMCMC$x)
  } else {
    selection.df[,2*k-1] <- as.data.frame(dens.indRepMCMC$x)
  }
  selection.df[,2*k] <- dens.indRepMCMC$y
  #colnames(selection.df[c(2*k-1,2*k)]) <- c(paste0(candidate_regions[k],"_selection_time"), paste0(candidate_regions[k],"_selection_density"))
}

### Write out indrep MCMC results as csv table, for convenient plotting
setwd("/Users/luqman/Desktop/Migration time demography")
colnames(selection.df) <- c("Dsyl_CEN_selection_time", "Dsyl_CEN_selection_density", "Dsyl_TCF1_selection_time", "Dsyl_TCF1_selection_density", "Dcar_TCF1_selection_time", "Dcar_TCF1_selection_density", "Dcar_FT_selection_time", "Dcar_FT_selection_density" )
write.table(selection.df, file = "selectionTimeEstimates_Free_newSS.txt", row.names = FALSE)
