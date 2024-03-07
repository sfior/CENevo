###### This script specifies the best-fitting theoretical distributions to the empirical data ######

### For reference, see: 
# https://r-forge.r-project.org/scm/viewvc.php/*checkout*/www/fitdistrplusE.pdf?revision=19&root=riskassessment&pathrev=21
# https://cran.r-project.org/web/packages/fitdistrplus/vignettes/paper2JSS.pdf
# http://stats.stackexchange.com/questions/132652/how-to-determine-which-distribution-fits-my-data-best


### First install and load packages
#install.packages("fitdistrplus")
#install.packages("logspline")
#install.packages("MASS")
#install.packages("qpcR")
library(fitdistrplus)
library(logspline)
library(MASS)
library(ggplot2)
library(qpcR)

### Import data (DEPTH PER SITE format)
# # This data is depreciated, but represents working example
# # test file with 1000 observations
# raw_data <- read.table("/Users/luqman/Documents/Hirzi/ETHPHD/Analysis/Coverage per site/test.txt", header = FALSE, sep = "\t", row.names = NULL)
# # full raw coverage file
# raw_data <- read.table("/Users/luqman/Documents/Hirzi/ETHPHD/Analysis/Coverage per site/FBcalls.SNPs.gdepth.txt", header = FALSE, sep = "\t", row.names = NULL)
# # filtered coverage file (min_cov = 10, max_cov = 95; as per according to Simone's settings)
# raw_data <- read.table("/Users/luqman/Documents/Hirzi/ETHPHD/Analysis/Coverage per site/FBcalls.SNPs.gdepth.filtered.txt", header = FALSE, sep = "\t", row.names = NULL)
# 
# ### For testing, you may want to subsample data (since there are >17 million rows)
# raw_data_subsampled <- raw_data[sample(nrow(raw_data), 100000), ]

### Import data (DEPTH COUNTS format)
species <- "Dcar"
# Set working directory
if (species == "Dsyl") {
  setwd('/Users/luqman/Documents/Hirzi/ETHPHD/Analysis/Coverage per site/Dsyl_angsd_output_4_Hirzi')
  barcodes<-c('ACAGTG','CAGATC','CTTGTA','AGTCAA','GCCAAT','GTGAAA')
} else if (species == "Dcar") {
  setwd('/Users/luqman/Documents/Hirzi/ETHPHD/Analysis/Coverage per site/Dcar_angsd_output_4_Hirzi')
  barcodes<-c('CGTACG','GTGGCC','GTTTCG','ATTCCT','GATCAG','TTAGGC')
}
## Note the barcode lists above have been arranged according to the population arrangement in thhe poolsim pipeline.
# Dsyl
#HIGH_ACAGTG, HIGH_CAGATC, HIGH_CTTGTA, LOW_AGTCAA, LOW_GCCAAT, LOW_GTGAAA
# Dcar
#HIGH_CGTACG, HIGH_GTGGCC, HIGH_GTTTCG, LOW_ATTCCT, LOW_GATCAG, LOW_TTAGGC

# Generate depth per site dataframe by making n length vectors (n is the count) of the depth value m.
raw_data_subsampled<-c()
for (b in barcodes){
  v<-scan(paste0(b,'_bamfiles.param1.qc.depthGlobal'))
  # Recall that here, depths above 1200 have been binned into the final (1200th) element, so we remove it here
  v<-v[-length(v)]
  x<-c()
  for (n in 0:1199){
    m<-n+1
    # Let's downsample by 3000. Recall that we have circa 343 million sites, so this will give us ~114,000 sites. This also ensures that the histogram is not truncated sharply/prematurely, as the final values in the depth counts (~1500 counts for depth 1199) when divided by 3000 will equate to 0.5, hence all higher depths will have integer counts of effectively zero.
    a<-rep(n,round(v[m]/3000,0))
    x<-c(x,a)
  }
  raw_data_subsampled <- qpcR:::cbind.na(raw_data_subsampled, x)
}
raw_data_subsampled<-raw_data_subsampled[,-1]
# Add population labels
colnames(raw_data_subsampled)<-barcodes
raw_data_subsampled <- as.data.frame(raw_data_subsampled)

### Columns 3-8 (depthPerSite format) and 1-6 (depthCounts format) represent the coverages for the 6 populations in dataframe. Let's test with one population/column first.
pop_data <- raw_data_subsampled$TTAGGC
pop_data_df <- as.data.frame(pop_data)
pop_data_df <- na.omit(pop_data_df)
pop_data_full <- pop_data_df[pop_data_df$pop_data >= 0, ]

# For initial visualisation purposes, it can be useful to remove right-hand outliers/extreme values
pop_data_trunc <- pop_data_df[pop_data_df$pop_data <= 200, ]

### Plot histogram and CDF plot of data
par(mfrow=c(1,1))
plotdist(pop_data_full, histo = TRUE, demp = TRUE)
plotdist(pop_data_trunc, histo = TRUE, demp = TRUE)

### Characterise and explore which theoretical distribution best fits data
#help(descdist)
# Assuming discrete
#descdist(pop_data_full, discrete = TRUE)
#descdist(pop_data_trunc, discrete = TRUE)
descdist(pop_data_full, discrete = TRUE, boot = 100)
# Assuming continous
#descdist(pop_data)
#descdist(pop_data_trunc)
#descdist(pop_data, boot = 100)

### Fit the distribution to the data
pop_data_norm <- fitdist(pop_data_trunc, "norm", method="mme")
#pop_data_pois <- fitdist(pop_data_trunc, "pois", method="mme")
pop_data_nbinom <- fitdist(pop_data_trunc, "nbinom", method="mme")
#plot(pop_data_norm)
summary(pop_data_norm)
#plot(pop_data_pois)
#summary(pop_data_pois)
#plot(pop_data_nbinom)
summary(pop_data_nbinom)

# Compare theoretial distributions with empirical.
par(mfrow=c(2,2))
plot.legend <- c("normal", "negative binomial")
denscomp(list(pop_data_norm, pop_data_nbinom))
qqcomp(list(pop_data_norm, pop_data_nbinom))
cdfcomp(list(pop_data_norm, pop_data_nbinom))
ppcomp(list(pop_data_norm, pop_data_nbinom))

# Estimate uncertainty in parameter estimates
pop_data_norm_param_uncertainties <- bootdist(pop_data_norm, niter = 100)
summary(pop_data_norm_param_uncertainties)
#plot(pop_data_norm_param_uncertainties)
pop_data_nbinom_param_uncertainties <- bootdist(pop_data_nbinom, niter = 100)
summary(pop_data_nbinom_param_uncertainties)
#plot(pop_data_nbinom_param_uncertainties)


######### RESULTS #########

# DSYL - for all populations, the negative binomial is the best-fitting distribution. 
# Optimised parametric bootstrapped medians for the negative binomial distribution paramaters.
barcodes<-c('ACAGTG','CAGATC','CTTGTA','AGTCAA','GCCAAT','GTGAAA')
mu <- c(42.192217, 33.497154, 40.334565, 45.408570, 39.031634, 43.251578)
size <- c(2.216471, 1.911589, 2.117969, 2.334359, 2.108416, 2.239261)        
# Optimised parametric bootstrapped medians for the normal distribution paramaters.
mean <- c(42.17227, 33.47498, 40.33264, 45.37401, 39.04726, 43.24986)
sd <- c(29.04633, 24.91928, 28.41397, 30.47678, 27.60674, 29.65197)      

# DCAR - for all populations (apart from GTTTCG, where the normal is a better fit), the negative binomial is the best-fitting distribution. 
# Optimised parametric bootstrapped medians for the negative binomial distribution paramaters.
barcodes<-c('CGTACG','GTGGCC','GTTTCG','ATTCCT','GATCAG','TTAGGC')
mu <- c(55.164384, 51.700325, 61.988267, 51.603372, 49.716216, 48.471899)
size <- c(3.659316, 3.241173, 3.612606, 3.188717, 3.110145, 3.017108)        
# Optimised parametric bootstrapped medians for the normal distribution paramaters.
mean <- c(55.18763, 51.70975, 61.97955, 51.58287, 49.71583, 48.45615)
sd <- c(29.78927, 29.61423, 33.53293, 29.74324, 29.05367, 28.76135)

