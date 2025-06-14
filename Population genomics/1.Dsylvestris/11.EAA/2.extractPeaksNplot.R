library(tidyverse)

setwd("./")

# Load the data
df <- read.table('LFMM_results_Elevation_K2.txt.gz', header = TRUE)
df <- df %>% separate(snpid, into = c("scaffold", "pos"), sep = "_pos")

# Convert the data to a data frame
results.lfmm <- as.data.frame(df)

###### Adds the FDR corrected pvalues to the table
# Ensure adjpvalue column is numeric and handle NAs
results.lfmm$adjpvalue <- as.numeric(as.vector(results.lfmm$adjpvalue))
nrow(results.lfmm)
# Remove rows with NA p-values before applying the adjustment
results.lfmm <- results.lfmm[!is.na(results.lfmm$adjpvalue), ]
nrow(results.lfmm)
# Apply FDR correction using the BH procedure
results.lfmm$FDRpvalue <- p.adjust(results.lfmm$adjpvalue, method = "BH")
gz_file <- gzfile('FMM_results_Elevation_K2_wFRDpval.txt.gz', 'w')
write.table(results.lfmm, gz_file, row.names = FALSE, quote = FALSE, sep = "\t")
close(gz_file)

# Filter the results based on the adjusted p-values
# Set FDR threshold
fdr.thres <- 0.05
results.lfmm.sign <- results.lfmm[results.lfmm$FDRpvalue <= fdr.thres, ]
write.table(results.lfmm.sign, 'FMM_results_Elevation_K2_wFRDpval_q0.05.txt', row.names = FALSE, quote = FALSE, sep = "\t")

###### plot scaffold4 with FDR values
# Subset the data
sc4 <- subset(results.lfmm, scaffold == 'scaffold4_size532381')
# Ensure FDRpvalue is numeric
sc4$FDRpvalue <- as.numeric(sc4$FDRpvalue)
nrow(sc4)
# Remove NA values if any (optional)
sc4 <- na.omit(sc4)
nrow(sc4)
# Transform p-values to log scale
sc4$logFDRpvalue <- -log10(sc4$FDRpvalue)
# Create the plot with log-transformed p-values
plot(sc4$pos, sc4$logFDRpvalue, cex = 0.5, main = "Log-transformed FDR p-values for Scaffold 4",
     xlab = "Position", ylab = "-log10(FDR p-value)")
# Add horizontal lines (transformed thresholds)
abline(h = -log10(0.05), col = "red", lwd = 1)
abline(h = -log10(0.1), col = "red", lwd = 1)
abline(h = -log10(0.2), col = "red", lwd = 1)
# Add orange points for DI peak
black_points <- subset(sc4, pos > 179265 & pos < 204502)
points(black_points$pos, black_points$logFDRpvalue, col = "orange", pch = 19, cex = 0.5)
# Add red points for CEN gene
red_points <- subset(sc4, pos > 191977 & pos < 193338)
points(red_points$pos, red_points$logFDRpvalue, col = "red", pch = 19, cex = 0.5)


##### Associations in DI peaks
# Read candidate peaks
results.lfmm$pos<-as.numeric(as.vector(results.lfmm$pos))
# Retain associations within peaks
DIpeaks<-read.table('candidate_regions_DI.txt',header=T)
results.lfmm.peaks<-c()
for (n in 1:nrow(DIpeaks)){
  scaff<-DIpeaks$scaffold[n]
  start<-DIpeaks$start[n]
  end<-DIpeaks$end[n]
  tmp<-subset(results.lfmm,scaffold==scaff & pos > start & pos < end)
  results.lfmm.peaks<-rbind(results.lfmm.peaks,tmp)
}
write.table(results.lfmm.peaks, 'FMM_results_Elevation_K2_wFRDpval_DIpeaks.txt', row.names = FALSE, quote = FALSE, sep = "\t")

# Extract the strongest association within each peak
results.lfmm.1top<-c()
for (n in 1:nrow(DIpeaks)){
  scaff<-DIpeaks$scaffold[n]
  tmp<-subset(results.lfmm.peaks,scaffold==scaff)
  # Find the index of the minimum value in the specified column
  min_index <- which.min(tmp[["FDRpvalue"]])
  # Subset the data frame to keep only the row with the minimum value
  top.snp <- tmp[min_index, ]
  # 
  results.lfmm.1top<-rbind(results.lfmm.1top,top.snp)
}
write.table(results.lfmm.1top, 'FMM_results_Elevation_K2_wFRDpval_DIpeaks_1top.txt', row.names = FALSE, quote = FALSE, sep = "\t")

# Plot signal for genome-wide LFMM with fdr 0.5 and add beta for top FDR SNP in each peak
data<-read.table('FMM_results_Elevation_K2_wFRDpval_q0.05.txt',header=T)
results.lfmm.1top<-read.table('FMM_results_Elevation_K2_wFRDpval_DIpeaks_1top.txt',header=T,sep='\t')
# scaffold 440 has 2 distinct regions, drop one (they have similar beta values)
results.lfmm.1top<-results.lfmm.1top[-4,]
# add number of candidates for plotting the values
results.lfmm.1top$N<-c(1:36)

hist(abs(data$beta)*1000,breaks = 200)
for (n in 1:nrow(results.lfmm.1top)){
  fdr<-results.lfmm.1top$FDRpvalue[n]
  beta<-results.lfmm.1top$beta[n]
  value_N <- results.lfmm.1top$N[n]  # Extract the value from column 'N'
  
  if (fdr < 0.05){col<-'red'}
  else {col<-'blue'}
  #abline(v = abs(beta)*1000, col = col, lwd = 2)
  # Plot small triangles below the x-axis
  x_pos <- abs(beta) * 1000
  points(x_pos, -1, col = col, pch = 17, cex = 1.5)  # pch = 17 is a triangle
  text(x_pos, 0, labels = value_N, col = 'black', pos = 3, cex = 0.8)
}  
sc4<-subset(results.lfmm.1top, scaffold=='scaffold4_size532381')
abline(v = abs(sc4$beta)*1000, col = 'red', lwd = 2)



