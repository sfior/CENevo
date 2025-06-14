### This script calculates, plots and outputs the averages (mean, median and mode) and credible intervals (highest density region and quantile) for the coalescent and speciation times inferred by Mathias.

# Load libraries
library("stringr")

# Set working directory
setwd("/Users/hl636/Documents/Hirzi/ETHPHD/Migration time demography")

# To convert number of generations to number of years, define the generation time
gen_time = 4

### Import data (my selection time estimates and Mathias' demography estimates)
# Speciation time estimates, Mathias
mathias_data <- read.table("exponential_growth.450k.posterior_distributions.txt", header = TRUE)
speciation_times <- mathias_data$T_merge_2_into_1
# Coalescence times estimates, Mathias
mathias_data2 <- read.table("coalescent_times.ddRAD_Dianthus.SI.exponential_growth.CEN_scheme_1column.txt", header = TRUE)
coalescence_times <- mathias_data2[,1]

### Transform data to appropriate format and scale
# Convert time to number of generations (in my demography analysis, Ne = 10,000; Mathias' demography analysis, Ne = 100,000), or if gen_time != 1, number of years
speciation_times <- speciation_times * 4 * 100000 * gen_time
coalescence_times <- coalescence_times * 4 * 100000 * gen_time
# Log transform
speciation_times <- log10(speciation_times)
coalescence_times <- log10(coalescence_times)
# Normalise
speciation_dens <- density(speciation_times, bw = "nrd", adjust = 1.5, kernel = "gaussian")
speciation_dens$y <- speciation_dens$y / max(speciation_dens$y)
coalescence_dens <- density(coalescence_times, bw = "SJ", adjust = 1, kernel = "gaussian")
coalescence_dens$y <- coalescence_dens$y / max(coalescence_dens$y)

# Format as dataframe
coal.df <- data.frame(coal_time=coalescence_dens$x, coal_dens=coalescence_dens$y, spec_time=speciation_dens$x, spec_dens=speciation_dens$y)

# Define credible interval parameters
# For reference, see: https://stats.stackexchange.com/questions/240749/how-to-find-95-credible-interval.
# Generally, quantiles will give you the interval containing probability mass concentrated around the median (the middle 100𝛼% of your distribution), while highest density region is a region around the modes of the distribution
credible_interval_1D <- 0.95

# First we initialise results dataframe, where we'll write the results to
CI_summary_df <- data.frame(param=character(), mean=character(), median=character(), mode=character(), lowerCI_HDR=character(), upperCI_HDR=character(), lowerCI_quantile=character(), upperCI_quantile=character(), stringsAsFactors=FALSE)

# Calculate distribution mean, median, mode and credible intervals (highest density region)
for (k in seq(1,3,2)) {
  print(k)
  param_times <- coal.df[,k]
  param_density <- coal.df[,k+1]
  
  ## Calculate CIs
  # Highest density region CI
  const <- sum(param_density)
  spxx <- sort(param_density, decreasing = TRUE) / const
  which(cumsum(spxx) >= credible_interval_1D)[1]
  HDR_CI <- spxx[which(cumsum(spxx) >= credible_interval_1D)[1]] * const
  # Get intersects (lower and upper)
  time_intersect_idx <- which.min(abs(param_density - HDR_CI))
  if(length(time_intersect_idx) < 2) {
    #time_intersect_idx <- sort(c(time_intersect_idx, which.min(abs(param_density[-time_intersect_idx] - HDR_CI))), decreasing = FALSE)
    time_intersect_idx <- sort(c(time_intersect_idx, which.min(abs(param_density[-c(time_intersect_idx-1, time_intersect_idx, time_intersect_idx + 1)] - HDR_CI))), decreasing = FALSE)
  }
  time_intersect_HDR <- param_times[time_intersect_idx]
  # Quantile interval
  cpxx <- cumsum(param_density) / sum(param_density)
  lower_bound_CI <- param_times[which(cpxx >= (0 + ((1 - credible_interval_1D) / 2)))[1]] 
  upper_bound_CI <- param_times[which(cpxx >= (1 - ((1 - credible_interval_1D) / 2)))[1]-1]
  time_intersect_quantile <- c(lower_bound_CI, upper_bound_CI)
  
  ## Calculate averages
  # Mean
  time_density_product <- param_times*param_density
  weighted_mean_dist <- sum(time_density_product)/sum(param_density)
  # Median
  cpxx <- cumsum(param_density) / sum(param_density)
  median_dist <- param_times[which(cpxx >= (0.5))[1]] 
  # Mode
  mode_dist <- param_times[which.max(param_density)]
  # Concatenate
  averages_dist <- c(weighted_mean_dist,median_dist, mode_dist)
  
  ## Write summary table
  candidate_label_long <- as.vector(str_split(colnames(coal.df[k]),"_",simplify = TRUE))
  candidate_label <- paste0(candidate_label_long[1],"_",candidate_label_long[2])
  data.row <- c(candidate_label,averages_dist,time_intersect_HDR, time_intersect_quantile)
  CI_summary_df[nrow(CI_summary_df) + 1, ] <- data.row
  
  # Plot
  titles <- c("Posterior distribution - speciation time (D.sylvestris - D.superbus)", "Posterior distribution - coalescence time (Dianthus sylvestris)")
  plot(param_times, param_density, type = "n", main=titles[1%%k + 1], ylab = " Posterior density", xlab = paste0("log10 years (assumed generation time: ", gen_time, " years per generation)"))
  polygon(param_times, param_density, col="darkgoldenrod3", border = NA)
  abline(v=time_intersect_HDR, lwd = 2, lty = 3, col="firebrick1")
  abline(v=time_intersect_quantile, lwd = 2, lty = 3, col="firebrick4")
  abline(v=averages_dist, lwd = 2, lty = 2, col=c("dodgerblue3", "darkorchid3", "olivedrab3"))
  #Add legend
  legend("topleft", inset = 0, c("mean", "median", "mode", "HDR CIs", "quantile CIs"), col=c("dodgerblue3", "darkorchid3", "olivedrab3", "firebrick1", "firebrick4"), cex=1, lty=c(2,2,2,3,3), lwd = 2, box.lty=0)
}

# Write out summary table
write.table(CI_summary_df, file = "coalTimeEstimates_averages_CIs.txt", row.names = FALSE, quote = FALSE, sep = "\t")
