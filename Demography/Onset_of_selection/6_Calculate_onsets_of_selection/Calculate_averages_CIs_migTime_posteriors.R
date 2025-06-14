### This script calculates, plots and outputs the averages (mean, median and mode) and credible intervals (highest density region and quantile) for the posterior probability distributions of candidate regions' time of onset of selection.
### It takes as input the selection time estimates produced via Selection_time_estimation_allRegions.R

# Load libraries
library("stringr")

# Set working directory
setwd("/Users/hl636/Documents/Hirzi/ETHPHD/Migration time demography")

# Import data (selection time estimates)
selection.df <- read.table("selectionTimeEstimates_Free_newSS.txt", header = TRUE) # simpleModel_11params_RECON2_migTime_FINAL_RUN2_FREE
candidate_regions <- c("Dsyl_CEN", "Dsyl_TCF1", "Dcar_TCF1", "Dcar_FT")

# To convert number of generations to number of years, define the generation time
gen_time = 4

# Define credible interval parameters
# For reference, see: https://stats.stackexchange.com/questions/240749/how-to-find-95-credible-interval.
# Generally, quantiles will give you the interval containing probability mass concentrated around the median (the middle 100𝛼% of your distribution), while highest density region is a region around the modes of the distribution
credible_interval_1D <- 0.95

# First we initialise results dataframe, where we'll write the results to
CI_summary_df <- data.frame(locus=character(), mean=character(), median=character(), mode=character(), lowerCI_HDR=character(), upperCI_HDR=character(), lowerCI_quantile=character(), upperCI_quantile=character(), stringsAsFactors=FALSE)

# Calculate distribution mean, median, mode and credible intervals (highest density region)
for (k in match(paste0(candidate_regions, "_selection_time"), colnames(selection.df))) {
  print(k)
  selection_times <- selection.df[,k]
  selection_density <- selection.df[,k+1]
  
  # Convert time to number of generations (in my demography analysis, Ne = 10,000; Mathias' demography analysis, Ne = 100,000), or if gen_time != 1, number of years
  selection_times <- (10^selection_times) * 4 * 10000 * gen_time
  # Log transform
  selection_times <- log10(selection_times)
  
  ## Calculate CIs
  # Highest density region CI
  const <- sum(selection_density)
  spxx <- sort(selection_density, decreasing = TRUE) / const
  which(cumsum(spxx) >= credible_interval_1D)[1]
  HDR_CI <- spxx[which(cumsum(spxx) >= credible_interval_1D)[1]] * const
  # Get intersects (lower and upper)
  time_intersect_idx <- which.min(abs(selection_density - HDR_CI))
  if(length(time_intersect_idx) < 2) {
    time_intersect_idx <- sort(c(time_intersect_idx, which.min(abs(selection_density[-time_intersect_idx] - HDR_CI))), decreasing = FALSE)
  }
  time_intersect_HDR <- selection_times[time_intersect_idx]
  # Quantile interval
  cpxx <- cumsum(selection_density) / sum(selection_density)
  lower_bound_CI <- selection_times[which(cpxx >= (0 + ((1 - credible_interval_1D) / 2)))[1]] 
  upper_bound_CI <- selection_times[which(cpxx >= (1 - ((1 - credible_interval_1D) / 2)))[1]-1]
  time_intersect_quantile <- c(lower_bound_CI, upper_bound_CI)
  
  ## Calculate averages
  # Mean
  time_density_product <- selection_times*selection_density
  weighted_mean_dist <- sum(time_density_product)/sum(selection_density)
  # Median
  cpxx <- cumsum(selection_density) / sum(selection_density)
  median_dist <- selection_times[which(cpxx >= (0.5))[1]] 
  # Mode
  mode_dist <- selection_times[which.max(selection_density)]
  # Concatenate
  averages_dist <- c(weighted_mean_dist,median_dist, mode_dist)
  
  ## Write summary table
  candidate_label_long <- as.vector(str_split(colnames(selection.df[k]),"_",simplify = TRUE))
  candidate_label <- paste0(candidate_label_long[1],"_",candidate_label_long[2])
  data.row <- c(candidate_label,averages_dist,time_intersect_HDR, time_intersect_quantile)
  CI_summary_df[nrow(CI_summary_df) + 1, ] <- data.row
  
  # Plot
  plot(selection_times, selection_density, type = "n", main=paste0("Posterior distribution for onset of selection: ", candidate_label), ylab = " Posterior density", xlab = paste0("log10 years (assumed generation time: ", gen_time, " years per generation)"), xlim = c(3,7.5))
  polygon(selection_times, selection_density, col="darkgoldenrod3", border = NA)
  abline(v=time_intersect_HDR, lwd = 2, lty = 3, col="firebrick1")
  abline(v=time_intersect_quantile, lwd = 2, lty = 3, col="firebrick4")
  abline(v=averages_dist, lwd = 2, lty = 2, col=c("dodgerblue3", "darkorchid3", "olivedrab3"))
  #Add legend
  legend("topright", inset = 0, c("mean", "median", "mode", "HDR CIs", "quantile CIs"), col=c("dodgerblue3", "darkorchid3", "olivedrab3", "firebrick1", "firebrick4"), cex=1, lty=c(2,2,2,3,3), lwd = 2, box.lty=0)
}

# Write out summary table
write.table(CI_summary_df, file = "selectionTimeEstimates_Free_newSS_averages_CIs.txt", row.names = FALSE, quote = FALSE, sep = "\t")
