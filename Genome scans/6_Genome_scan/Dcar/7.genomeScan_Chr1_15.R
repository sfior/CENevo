############# This script performs a genome scan based on demographic parameter posteriors of sliding windows in a chromosome. #############
## I.e. it identifies outlier windows based on inferred demographic parameters that deviate from neutrality
## It takes as the input file the output of ABCestimator
## This version performs a genome scan even if (the regression in ABCestimator) fails for some windows (e.g. no data, 0 segsites). It does this by removing these failed windows from the scan.

# Read arguments from command line
args <- commandArgs(trailingOnly = TRUE)

# Import libraries
library(MASS)
library(stringr)
library(gtools)

############# Define input variables #############

# Don't forget to adjust the priors between lines 142-143 accordingly!

# Input directory (directory containing ABC estimation results)
#input_dir <- "/Users/luqman/Desktop/ABC GS Revisions/"
#input_dir <- "/cluster/scratch/lhirzi/ABCEstimator"
input_dir <- "/cluster/project/gdc/shared/p219/ABCtoolbox_run/Estimation_results"
# Input directory prefix
prefix_dir <- "Estimation_results_simpleModel_11params_RECON2_newSS_minDP10maxDP200_Dcar_3RunsCombined_PLS20_Joint_retSims10000_gridPoints33_FINAL"
# Species (Dsyl or Dcar)
species <- "Dcar"
# Chromosome
#Chr <- 1
# Number of parameters in demographic model
nparams <- 11
# Define whether joint posteriors were calculated via grid search "grid" or MCMC sampling "MCMC"
joint.type = "grid"
# Parameter names
param_labels <- c("mIC", "mHL", "mLH", "NH", "NL", "NI_L1", "NI_L2", "NI_L3", "NI_H1", "NI_H2", "NI_H3")
# Credible interval to define significance
credible_interval <- 0.999
#credible_interval_1D <- round(credible_interval^0.5,3) # remember, for equivalence we take the square root of the 2D CI for the 1D CI
credible_interval_1D <- 0.995

if (species == "Dsyl") {
  # Marginal neutral estimates taken from Master_plotter_simpleModel_11params_XXX.R
  neutral_estimate <- c(3.33333, 2.57576, 2.51515, 5.82828, 5.78789, 3.57576, 3.22222, 3.0202, 2.56566, 2.81818, 2.21212) # RECON2 newSS minDP10maxDP200 3RunsCombined NEW
  # Joint neutral estimates taken from Master_plotter_simpleModel_11params_XXX.R (marginalNeutral) or from Combine_posteriors_2DJoint.R (jointNeutral). The former is preferred because 1) it gives a better signal and 2) it has been validated through simulations with pseudo-observed data.
  #neutral_estimate_2D <- c(2.63636, 0.818182) # RECON2_newSS_minDP10maxDP200_Dsyl_3RunsCombinedjointNeutral Grid
  neutral_estimate_2D <- c(2.57576, 2.51515) # RECON2_newSS_minDP10maxDP200_Dsyl_3RunsCombined marginalNeutral 
  # Number of windows per chromosome
  num_loci_list <- c(6665, 4792, 7620, 5568, 5750, 6412, 4855, 4950, 5405, 5430, 4249, 4245, 7437, 7690, 4963, 135246) #RECON2 newSS 1kb overlapping windows zeroSegSitesRemoved  # List of scaffolds
  Chr_scaffold_list <- c("Chr1_allScaffolds_FullSS_start1_end7187273_windowsize5000_M2_Q0_1kbOverlapWindows.zeroSegSitesRemoved.scaffoldlist", "Chr2_allScaffolds_FullSS_start1_end5149401_windowsize5000_M2_Q0_1kbOverlapWindows.zeroSegSitesRemoved.scaffoldlist", "Chr3_allScaffolds_FullSS_start1_end8134089_windowsize5000_M2_Q0_1kbOverlapWindows.zeroSegSitesRemoved.scaffoldlist", "Chr4_allScaffolds_FullSS_start1_end5958733_windowsize5000_M2_Q0_1kbOverlapWindows.zeroSegSitesRemoved.scaffoldlist", "Chr5_allScaffolds_FullSS_start1_end6199402_windowsize5000_M2_Q0_1kbOverlapWindows.zeroSegSitesRemoved.scaffoldlist", "Chr6_allScaffolds_FullSS_start1_end6750011_windowsize5000_M2_Q0_1kbOverlapWindows.zeroSegSitesRemoved.scaffoldlist", "Chr7_allScaffolds_FullSS_start1_end5209322_windowsize5000_M2_Q0_1kbOverlapWindows.zeroSegSitesRemoved.scaffoldlist", "Chr8_allScaffolds_FullSS_start1_end5355932_windowsize5000_M2_Q0_1kbOverlapWindows.zeroSegSitesRemoved.scaffoldlist", "Chr9_allScaffolds_FullSS_start1_end5836137_windowsize5000_M2_Q0_1kbOverlapWindows.zeroSegSitesRemoved.scaffoldlist", "Chr10_allScaffolds_FullSS_start1_end5817725_windowsize5000_M2_Q0_1kbOverlapWindows.zeroSegSitesRemoved.scaffoldlist", "Chr11_allScaffolds_FullSS_start1_end4529089_windowsize5000_M2_Q0_1kbOverlapWindows.zeroSegSitesRemoved.scaffoldlist", "Chr12_allScaffolds_FullSS_start1_end4630051_windowsize5000_M2_Q0_1kbOverlapWindows.zeroSegSitesRemoved.scaffoldlist", "Chr13_allScaffolds_FullSS_start1_end7920842_windowsize5000_M2_Q0_1kbOverlapWindows.zeroSegSitesRemoved.scaffoldlist", "Chr14_allScaffolds_FullSS_start1_end8187860_windowsize5000_M2_Q0_1kbOverlapWindows.zeroSegSitesRemoved.scaffoldlist", "Chr15_allScaffolds_FullSS_start1_end5372557_windowsize5000_M2_Q0_1kbOverlapWindows.zeroSegSitesRemoved.scaffoldlist", "Chr16_allScaffolds_FullSS_start1_end9999999_windowsize5000_M2_Q0_1kbOverlapWindows.zeroSegSitesRemoved.scaffoldlist")  #RECON2 newSS 1kb overlapping windows zeroSegSitesRemoved
  # Import data (linkage group position and chromosome number)
  setwd("/cluster/project/gdc/people/lhirzi/Assembly/Dsylvestris/anchored_scaffolds")
  #setwd("/Users/luqman/Desktop/GenomeScan")
  scaffold_pos_raw<-read.delim("Dsyl_all_anchored_scaffolds_v2.txt", sep = "\t", header = TRUE)
  # Output prefix
  output_dir <- "/cluster/project/gdc/people/lhirzi/GenomeScan/Dsyl/GenomeScan_resultsTables"
  #output_dir <- "/Users/luqman/Desktop"
} else if (species == "Dcar") {
  # Neutral estimates taken from Master_plotter_simpleModel_11params_RECON2.R
  neutral_estimate <- c(2.87879, 2.39393, 2.51515, 5.50505, 5.50508, 3.17172, 3.22222, 3.22222, 2.11112, 2.11111, 2.36364) # RECON2 newSS minDP10maxDP200 3RunsCombined NEW
  neutral_estimate_2D <- c(2.39393, 2.51515) # RECON2_newSS_minDP10maxDP200_Dsyl_3RunsCombined - taking the marginals
  # Number of windows per chromosome
  num_loci_list <- c(5187, 5296, 5208, 4927, 4448, 3927, 3782, 6826, 7560, 4884, 6604, 5396, 5546, 4789, 4610, 85904) #RECON2 newSS 1kb overlapping windows zeroSegSitesRemoved  # List of scaffolds
  # List of scaffolds
  Chr_scaffold_list <- c("Chr1_allScaffolds_FullSS_start1_end5522201_windowsize5000_M2_Q0_1kbOverlapWindows.zeroSegSitesRemoved.scaffoldlist", "Chr2_allScaffolds_FullSS_start1_end5623217_windowsize5000_M2_Q0_1kbOverlapWindows.zeroSegSitesRemoved.scaffoldlist", "Chr3_allScaffolds_FullSS_start1_end5547344_windowsize5000_M2_Q0_1kbOverlapWindows.zeroSegSitesRemoved.scaffoldlist", "Chr4_allScaffolds_FullSS_start1_end5248954_windowsize5000_M2_Q0_1kbOverlapWindows.zeroSegSitesRemoved.scaffoldlist", "Chr5_allScaffolds_FullSS_start1_end4732141_windowsize5000_M2_Q0_1kbOverlapWindows.zeroSegSitesRemoved.scaffoldlist", "Chr6_allScaffolds_FullSS_start1_end4217913_windowsize5000_M2_Q0_1kbOverlapWindows.zeroSegSitesRemoved.scaffoldlist", "Chr7_allScaffolds_FullSS_start1_end4006302_windowsize5000_M2_Q0_1kbOverlapWindows.zeroSegSitesRemoved.scaffoldlist", "Chr8_allScaffolds_FullSS_start1_end7220078_windowsize5000_M2_Q0_1kbOverlapWindows.zeroSegSitesRemoved.scaffoldlist", "Chr9_allScaffolds_FullSS_start1_end7995774_windowsize5000_M2_Q0_1kbOverlapWindows.zeroSegSitesRemoved.scaffoldlist", "Chr10_allScaffolds_FullSS_start1_end5278334_windowsize5000_M2_Q0_1kbOverlapWindows.zeroSegSitesRemoved.scaffoldlist", "Chr11_allScaffolds_FullSS_start1_end6967360_windowsize5000_M2_Q0_1kbOverlapWindows.zeroSegSitesRemoved.scaffoldlist", "Chr12_allScaffolds_FullSS_start1_end5814850_windowsize5000_M2_Q0_1kbOverlapWindows.zeroSegSitesRemoved.scaffoldlist", "Chr13_allScaffolds_FullSS_start1_end5930082_windowsize5000_M2_Q0_1kbOverlapWindows.zeroSegSitesRemoved.scaffoldlist", "Chr14_allScaffolds_FullSS_start1_end5148409_windowsize5000_M2_Q0_1kbOverlapWindows.zeroSegSitesRemoved.scaffoldlist", "Chr15_allScaffolds_FullSS_start1_end4927259_windowsize5000_M2_Q0_1kbOverlapWindows.zeroSegSitesRemoved.scaffoldlist", "Chr16_allScaffolds_FullSS_start1_end9999999_windowsize5000_M2_Q0_1kbOverlapWindows.zeroSegSitesRemoved.scaffoldlist")  #RECON2 newSS 1kb overlapping windows zeroSegSitesRemoved
  # Import data (linkage group position and chromosome number)
  setwd("/cluster/project/gdc/people/lhirzi/Assembly/Dcarthusian/anchored_scaffolds")
  #setwd("/Users/luqman/Desktop/GenomeScan")
  scaffold_pos_raw<-read.delim("Dcart_all_anchored_scaffolds_v2.txt", sep = "\t", header = TRUE)
  # Output prefix
  output_dir <- "/cluster/project/gdc/people/lhirzi/GenomeScan/Dcar/GenomeScan_resultsTables"
  #output_dir <- "/Users/luqman/Desktop"
}

############# Run Genome Scan for nChr chromosomes #############

# Assign argument to chromosome number (to parallelise across CPUs)
Chr <- as.integer(args[1])

# Define number of windows or loci in chromosome
num_loci <- num_loci_list[Chr]

# Set working directory to point to marginalPosteriorDensities results
setwd(paste0(input_dir,"/",prefix_dir,"/",prefix_dir,"_Chr", Chr,"/marginalPosteriorDensities"))

# In case regression in ABCestimator fails for some loci (which can happen), we'll iterate over all succesfully estimated loci rather than all observations.
ABC_output_files <- mixedsort(list.files())
idx_loci_passed <- grep("MarginalPosteriorDensities", ABC_output_files)
# Check that indexes are sorted as expected
sort_result <- !is.unsorted(idx_loci_passed)
print(paste0("Loci are sorted by observation: ", sort_result))
# Define number of windows or loci in chromosome
#num_loci <- num_loci_list[Chr]
num_loci_passed <- length(idx_loci_passed)

############# Calculate posterior density peaks and marginal (1D) confidence intervals #############
posterior_peaks_matrix <- matrix(0, ncol=nparams, nrow=num_loci_passed)
#CI_matrix <- matrix(0, ncol=nparams, nrow=num_loci_passed)
significance_HDR_matrix <- matrix(0, ncol=nparams, nrow=num_loci_passed)
significance_QI_matrix <- matrix(0, ncol=nparams, nrow=num_loci_passed)
neutral_CI_matrix <- matrix(0, ncol=nparams, nrow=num_loci_passed)
#neutral_density_matrix <- matrix(0, ncol=nparams, nrow=num_loci_passed)
counter = 1
for (i in ABC_output_files[idx_loci_passed]) {
  for(p in 1:nparams){
    ABC_GLM<-read.delim(paste(i, sep = "")) 
    # Find peaks of posterior
    posterior_peak <- ABC_GLM[,2*p][match(max(ABC_GLM[,2*p+1]),ABC_GLM[,2*p+1])] # Probably want to use a better estimator than the max?
    # Output results to matrix
    posterior_peaks_matrix[counter,p] <- posterior_peak
    # Calculate 1D marginal credible interval (highest density region and quantile interval)
    # Highest density region CI
    const <- sum(ABC_GLM[,2*p+1])
    spxx <- sort(ABC_GLM[,2*p+1], decreasing = TRUE) / const
    HDR_CI <- spxx[which(cumsum(spxx) >= credible_interval_1D)[1]] * const
    # % credible interval at neutral estimate (as a continuous measure of significance)
    neutral_density_index <- which.min(abs(ABC_GLM[,2*p] - neutral_estimate[p]))
    neutral_density <- ABC_GLM[,2*p+1][neutral_density_index]
    neutral_index <- match((neutral_density / const), spxx)
    neutral_CI <- cumsum(spxx)[neutral_index]   # this variable indicates at which HDR credible interval the neutral estimate lies (in fraction)
    neutral_CI_matrix[counter,p] <- neutral_CI
    # Quantile interval
    cpxx <- cumsum(ABC_GLM[,2*p+1]) / sum(ABC_GLM[,2*p+1])
    lower_bound_CI <- ABC_GLM[,2*p][which(cpxx >= (0 + ((1 - credible_interval_1D) / 2)))[1]] 
    upper_bound_CI <- ABC_GLM[,2*p][which(cpxx >= (1 - ((1 - credible_interval_1D) / 2)))[1]-1]
    # Calculate density at neutral estimate
    neutral_closestObs <- which.min(abs(ABC_GLM[[2*p]] - neutral_estimate[p]))
    neutral_est_density <- ABC_GLM[[2*p+1]][neutral_closestObs]
    # Calculate whether neutral point lies within 1D credible interval
    neutrality_inf_HDR <- neutral_est_density > HDR_CI
    neutrality_inf_QI <- (upper_bound_CI > neutral_estimate[p]) & (neutral_estimate[p] > lower_bound_CI)
    # Output results to matrix
    #CI_matrix[counter,p] <- HDR_CI
    #neutral_density_matrix[counter,p] <- neutral_est_density
    significance_HDR_matrix[counter,p] <- neutrality_inf_HDR
    significance_QI_matrix[counter,p] <- neutrality_inf_QI
  }
  counter = counter + 1
}

############# Calculate 2D joint posterior confidence intervals #############
# It may be more relevant/powerful to draw credible intervals on the 2D joint posterior. We do this by performing a 2D KDE.
# Define run parameters
density_points = 33 # number of points per parameter.
prior_migHL <- c(-3,3) # RECON2
prior_migLH <- c(-3,3) # RECON2
# Initiate empty vectors
significance_KDE2D_vector <- vector()
significance_KDE2D_neutralCI_vector <- vector()
asymmetry_KDE2D_vector <- vector()
counter = 1

# Change working directory to point to jointPosterior results
setwd(paste0(input_dir,"/",prefix_dir,"/",prefix_dir,"_Chr", Chr,"/jointPosteriors"))

# In case regression in ABCestimator fails for some loci (which can happen), we'll iterate over all succesfully estimated loci rather than all observations.
ABC_output_files_joint <- mixedsort(list.files())
if (joint.type == "grid") {
  idx_loci_passed_joint <- grep("jointPosterior", ABC_output_files_joint)
} else if (joint.type == "MCMC") {
  idx_loci_passed_joint <- grep("jointPosteriorSamples", ABC_output_files_joint)
}
# Check that indexes are sorted as expected
sort_result_joint <- !is.unsorted(idx_loci_passed_joint)
print(paste0("Loci are sorted by observation: ", sort_result_joint))
# Define number of windows or loci in chromosome
#num_loci <- num_loci_list[Chr]
num_loci_passed_joint <- length(idx_loci_passed_joint)

for (i in ABC_output_files_joint[idx_loci_passed_joint]) {
  ABC_GLM<-read.delim(i, sep = "")
  # Specify whether joint (mHL and mLH) parameters were calculated via grid search or MCMC sampling
  if (joint.type == "grid") {
    dens <- list()
    # This (x-axis) is refers to log_m_highcontinent_lowcontinent
    dens[[1]] <- seq(prior_migHL[1],prior_migHL[2], length.out = density_points)
    # This (y-axis) is refers to log_m_lowcontinent_highcontinent
    dens[[2]] <- seq(prior_migLH[1],prior_migLH[2], length.out = density_points)
    # We transform the density vector into matrix.
    dens[[3]] <- matrix(ABC_GLM$density,nrow=density_points, ncol = density_points)
    names(dens) <- c("x", "y", "z")
  } else if (joint.type == "MCMC") {
    dens <- kde2d(ABC_GLM[,2], ABC_GLM[,3], n = density_points)
  }
  # Calculate and output asymmetry
  if (joint.type == "grid") {
    # Recall that there are points above the diagonal, below the diagonal and on the diagonal. To reduce degrees of freedom, we first account for the points on the diagonal (i.e. where mLH=Mhl).
    symmetry <- sum((ABC_GLM$log_m_highcontinent_lowcontinent == ABC_GLM$log_m_lowcontinent_highcontinent)*ABC_GLM$density)
    # Then the integral of the area above the diagonal equals one minus the integral of the area below the diagonal. So we can just output one variable.
    asymmetry <- sum((ABC_GLM$log_m_highcontinent_lowcontinent < ABC_GLM$log_m_lowcontinent_highcontinent)*ABC_GLM$density) / (sum(ABC_GLM$density) - symmetry)
  } else if (joint.type == "MCMC") {
    # Recall that there are points above the diagonal, below the diagonal and on the diagonal. To reduce degrees of freedom, we first account for the points on the diagonal (i.e. where mLH=Mhl).
    symmetry <- sum(ABC_GLM$log_m_highcontinent_lowcontinent == ABC_GLM$log_m_lowcontinent_highcontinent)
    # Then the integral of the area above the diagonal equals one minus the integral of the area below the diagonal. So we can just output one variable.
    asymmetry <- sum(ABC_GLM$log_m_highcontinent_lowcontinent < ABC_GLM$log_m_lowcontinent_highcontinent) / (length(ABC_GLM$log_m_highcontinent_lowcontinent) - symmetry)
  }
  # Append asymmetry statistic to results vector
  asymmetry_KDE2D_vector[counter] <- asymmetry
  # Output intersection point of neutral estimate and 2D joint distribution
  neutral_contour_x <- which.min(abs(dens$x - neutral_estimate_2D[1]))
  neutral_contour_y <- which.min(abs(dens$y - neutral_estimate_2D[2]))
  neutral_contour_xy <- dens$z[neutral_contour_x, neutral_contour_y]
  # Calculate credible interval (here: highest density interval)
  const2d <- sum(dens$z)
  spxx2d <- sort(dens$z, decreasing = TRUE) / const2d
  crit <- spxx2d[which(cumsum(spxx2d) >= credible_interval)[1]] * const2d
  # Append significance to results vector
  neutrality_inf_KDE2D <- neutral_contour_xy > crit
  significance_KDE2D_vector[counter] <- neutrality_inf_KDE2D
  # Additionally, we indicate at which % credible interval the neutral estimate lies (i.e. the area above the horizontal threshold defined by the neutral estimate)
  neutral_index2d <- match((neutral_contour_xy / const2d), spxx2d)
  neutral_CI_2d <- cumsum(spxx2d)[neutral_index2d]   # this variable indicates at which HDR credible interval the neutral estimate lies (in fraction)
  significance_KDE2D_neutralCI_vector[counter] <- neutral_CI_2d
  counter = counter + 1
}
significance_KDE2D_vector <- as.numeric(significance_KDE2D_vector)

############# For ploting #############
# We want to colour the points based on the credible interval overlap. We will do this by subsetting on matrix significance_matrix, which is currently defined as 1 = TRUE and 0 = FALSE). Since R indexes from 1 (not from 0), we add 1 to all elements of the matrix.
significance_HDR_matrix <- significance_HDR_matrix + 1
significance_QI_matrix <- significance_QI_matrix + 1
significance_KDE2D_vector <- significance_KDE2D_vector + 1

# To color points (windows) by their significance values (% credible interval where the neutral estimate lies), we first rank the variable of interest for colour assignment
significance_KDE2D_neutralCI_df <- as.data.frame(significance_KDE2D_neutralCI_vector)
significance_KDE2D_neutralCI_df$order = findInterval(significance_KDE2D_neutralCI_df$significance_KDE2D_neutralCI_vector, sort(significance_KDE2D_neutralCI_df$significance_KDE2D_neutralCI_vector))

# So that the lower values occupy more space, we transform the data as such: -log10(1-CI)
neutral_logCI_matrix <- -log10(1 - neutral_CI_matrix)
significance_KDE2D_neutralCI_df[,3] <- -log10(1 - significance_KDE2D_neutralCI_df[,1])

############# Make and output genomeScan table #############

# Read in scaffold list information
if (species == "Dsyl") {
  #setwd("/Users/luqman/Desktop/ABCtoolbox/")
  setwd(paste0("/cluster/project/gdc/people/lhirzi/GenomeScan/Dsyl/Chr_", Chr))
} else if (species == "Dcar") {
  setwd(paste0("/cluster/project/gdc/people/lhirzi/GenomeScan/Dcar/Chr_", Chr))
  #setwd("/Users/luqman/Desktop/ABCtoolbox/")
}
window_pos<-read.delim(Chr_scaffold_list[Chr], sep = "S", header = FALSE)

# Remove scaffolds if no posterior was calculated (i.e. if regression failed during ABCestimator)
first_passed_obs.raw <- ABC_output_files[idx_loci_passed[1]]
first_passed_obs.temp <- str_split_fixed(first_passed_obs.raw, "Obs", 2)
first_passed_obs <- str_split_fixed(first_passed_obs.temp[2], ".txt", 2)
first_passed_obs <- as.numeric(first_passed_obs[1])
idx_loci_passed_MarginalPosteriorDensities <- idx_loci_passed - idx_loci_passed[1] + (1 + first_passed_obs)
first_passed_obs_joint.raw <- ABC_output_files_joint[idx_loci_passed_joint[1]]
first_passed_obs_joint.temp <- str_split_fixed(first_passed_obs_joint.raw, "Obs", 2)
first_passed_obs_joint <- str_split_fixed(first_passed_obs_joint.temp[2], ".txt", 2)
first_passed_obs_joint <- as.numeric(first_passed_obs_joint[1])
idx_loci_passed_jointPosteriors <- idx_loci_passed_joint - idx_loci_passed_joint[1] + (1 + first_passed_obs_joint)
# Merge (find union) of idx_loci_passed (marginal and joint)
idx_loci_passed_combined <- intersect(idx_loci_passed_MarginalPosteriorDensities,idx_loci_passed_jointPosteriors)
window_pos <- window_pos[idx_loci_passed_combined,]

# Select relevant data
window_pos<-as.vector(window_pos[,4])
# Remove leading and trailing (useless) characters
scaffold_list.temp <- lapply(window_pos, function(x) substring(x, 2, nchar(x)-4))
scaffold_list <- as.data.frame(t(as.data.frame(scaffold_list.temp)))
rownames(scaffold_list)<-c(1:length(scaffold_list[,1]))
# Split into scaffold name and window range
scaffold_df <- as.data.frame(str_split_fixed(scaffold_list$V1, "_window_", 2))
scaffold_df.windows <- as.data.frame(str_split_fixed(scaffold_df$V2, "_", 2)) 

# Define dataframe
# Append scaffold name
genomeScan_table <- as.data.frame(scaffold_df$V1)
# Append scaffold window ranges (start and end)
genomeScan_table[,2] <- scaffold_df.windows$V1
genomeScan_table[,3] <- scaffold_df.windows$V2
# Append CI values (1D and 2D)
for (p in 2:3) {
  genomeScan_table[,p+2] <- neutral_logCI_matrix[,p]
}
genomeScan_table[,6] <- significance_KDE2D_neutralCI_df[,3]
# Append binary significance states based on 95% CI - 2:not significant, 1: significant (1D and 2D)
for (p in 2:3) {
  genomeScan_table[,p+5] <- significance_HDR_matrix[,p]
}
genomeScan_table[,9] <- significance_KDE2D_vector
# Append migration rate posterior modes 
for (p in 2:3) {
  genomeScan_table[,p+8] <- posterior_peaks_matrix[,p]
}
# Append asymmetry stastistic
genomeScan_table[,12] <- asymmetry_KDE2D_vector

# Add linkage group position and chromosome number
genomeScan_table[,13] <- match(genomeScan_table$scaffold, scaffold_pos_raw$scaffold)
genomeScan_table[,14] <- scaffold_pos_raw$postion[genomeScan_table[,13]]
genomeScan_table[,15] <- scaffold_pos_raw$LG[genomeScan_table[,13]]

# Add column headers
colnames(genomeScan_table) <- c("scaffold", "window_start", "window_end", "-log(1-p_migHL)", "-log(1-p_migLH)", "-log(1-p_joint_migHL_migLH)", paste0("significance_",credible_interval_1D,"_migHL"), paste0("significance_",credible_interval_1D,"_migLH"), paste0("significance_",credible_interval,"_joint_migHL_migLH"), "mode_migHL", "mode_migLH", "asymmetry", "pos_index", "Chr_position", "Chr")
# Reorder dataframe columns
genomeScan_table <- subset(genomeScan_table, select =  c("scaffold", "Chr", "Chr_position", "window_start", "window_end", "-log(1-p_migHL)", "-log(1-p_migLH)", "-log(1-p_joint_migHL_migLH)", paste0("significance_",credible_interval_1D,"_migHL"), paste0("significance_",credible_interval_1D,"_migLH"), paste0("significance_",credible_interval,"_joint_migHL_migLH"), "mode_migHL", "mode_migLH", "asymmetry") )
# Redefine column types. To convert factor to numeric, need to convert to character first (because factors are stored internally as integers with a table to give the factor level labels)
#sapply(genomeScan_table, class)
genomeScan_table$window_start <- as.numeric(as.character(genomeScan_table$window_start))
genomeScan_table$window_end <- as.numeric(as.character(genomeScan_table$window_end))

# Sort by chromosome position (first) and windows (second)
genomeScan_table <- genomeScan_table[with(genomeScan_table, order(Chr_position, window_start)), ]

# Write out table
setwd(output_dir)
write.table(genomeScan_table, file = paste0("genomeScan_table.Chr",Chr,".txt"), row.names = FALSE)
