
##### This script plots the selection time estimates for multiple candidate regions, in a ridge line plot, together with the estimate of speciation and coalescence times #####
## Note 1: It takes 3 files as inputs: 1) my selection time estimates (produced via Selection_time_estimation_allRegions.R) , 2) Mathias' speciation time estimates, and 3) Mathias' coalescence time estimates
## Note 2: Dsyl_CEN, Dsyl_TCF1, Dcar_TCF1, Dcar_FT posteriors were produced from 5, 2, 3, 1 windows respectively (take this into consideration)

### Define input variables
# Set working directory
setwd("/Users/luqman/Desktop/Migration time demography")
# Define which candidate regions to plot (any combination or number of: "Dsyl_CEN", "Dsyl_TCF1", "Dcar_TCF1", "Dcar_FT")
candidate_regions <- c("Dsyl_CEN", "Dsyl_TCF1", "Dcar_TCF1", "Dcar_FT")
#candidate_regions <- c("Dsyl_CEN", "Dsyl_TCF1", "Dcar_TCF1")
#candidate_regions <- c("Dsyl_TCF1", "Dcar_TCF1")
# Number of candidate regions assessed
num_regions = length(candidate_regions)
# To convert number of generations to number of years, define the generation time
gen_time = 4

### Import data (my selection time estimates and Mathias' demography estimates)
# Selection time estimates, Hirzi, independent replicates MCMC
selection.df <- read.table("selectionTimeEstimates_Free_newSS.txt", header = TRUE) # simpleModel_11params_RECON2_migTime_FINAL_RUN2_FREE
#selection.df <- read.table("selectionTimeEstimates_Fixed.txt", header = TRUE) # simpleModel_11params_RECON2_migTime_FINAL_RUN2_FIXED
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

### Plot posteriors
par(mfrow=c(1,1), mar=c(5.1,6.1,4.1,2.1))
# Plot background (axes, labels, title)
plot(speciation_dens, type='n', main="Speciation and selection times", yaxt = "n", ylab = NA , xlab = paste0("log10 years (assumed generation time: ", gen_time, " years per generation)"), xlim = c(4.5,7.5), ylim = c(0.175,0.9*(num_regions+1)));
axis(2, at = (seq(0.8,num_regions*0.8,0.8)), labels = NA, tick = TRUE)
text(y = (seq(0.8,num_regions*0.8,0.8)), par("usr")[1], labels = candidate_regions, srt = 0, pos = 2, xpd = TRUE)

# Plot selection time estimates for all candidate regions
counter = num_regions
for (k in rev(match(paste0(candidate_regions, "_selection_time"), colnames(selection.df)))) {
  selection_times <- selection.df[,k]
  selection_density <- selection.df[,k+1]
  # Convert time to number of generations (in my demography analysis, Ne = 10,000; Mathias' demography analysis, Ne = 100,000), or if gen_time != 1, number of years
  selection_times <- (10^selection_times) * 4 * 10000 * gen_time
  # Log transform
  selection_times <- log10(selection_times)
  # Plot
  #lines(selection_times, selection_density+counter-counter*(0.2), type='l', lwd = 1, col = rgb(counter*(1/(num_regions+1)), counter*(1/(num_regions+1)), 1, 1), lty = 1)
  line(selection_times, selection_density)
  polygon(selection_times, selection_density+counter-counter*(0.2), col=rgb(counter*(1/(num_regions+1)), counter*(1/(num_regions+1)), 1, 1), border = NA)
  # Update counter
  counter = counter - 1
}
# Plot speciation and coalescence time estimates 
line(speciation_dens)
polygon(speciation_dens,col=rgb(0.8, 0, 0, 0.75), border=NA)
line(coalescence_dens)
polygon(coalescence_dens,col=rgb(0, 0.8, 0, 0.75), border=NA)
# Plot time of radiation (taking value of 1.3-2 mya (Valente et al 2010))
# Assuming this coincides with 95% confidence interval (this is likely not true; the confidence intervals are much larger however not quantified (?) in the paper))
# x <- seq(0, 10, length=200)
# y <- dnorm(x, mean=1.65, sd=0.35)
# y <- y / max(y)
# x <- log10(x*1000000)
# Assuming larger confidence interval centered around 1.65 mya (center of 1.3-2 mya range)
x <- seq(0, 10, length=200)
y <- dnorm(x, mean=log10(1650000), sd=0.3)
y <- y / max(y)
line(x, y)  
polygon(x,y,col=rgb(0.95, 0.95, 0, 0.75), border=NA)
abline(v = log10(1650000), lty = 3, lwd = 2, col = "grey20")
#text(log10(2000000)*1.03,0.875*(num_regions+1),"Eurasian radiation event",srt=0.2,pos=3, col = "grey20", cex = 1)
#Add legend
legend("topright", inset = 0, c("onset of selection", "speciation time (D.sylvestris - D.superbus)", "coalescence time (Dianthus sylvestris)", "Eurasian radiation event"),
       fill=c("blue", "red3", "green3", "yellow"), cex=0.9, box.lty=0)

