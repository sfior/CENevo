
############# This script draws Manhatten plot for genome scans, for multiple (all) chromosomes plots #############

# Note 1: The input is a multiple (concatenated) chromomome genome scan table, in text format.
# Note 2: in the genome scan input table, a significance state of "1" implies significance and of "2" implies non-significant.
# Note 3: Before running, make sure to adjust the variables defined under plotting parameters at the top of the script.
# Note 4: This script additionally colors points above a defined threshold according to the symmetry of the posterior distribution.

# Load libraries
library(RColorBrewer)
library(DescTools)

# Define directory containing genome scan table files
#setwd("/Users/luqman/Desktop/GenomeScan/GenomeScan results/GenomeScan_resultsTables_simpleModel_11params_RECON2_newSS_minDP10maxDP200_Dcar_3RunsCombined_PLS20_Joint_retSims10000_gridPoints33_FINAL")
setwd("/Users/luqman/Desktop/GenomeScan/GenomeScan results/GenomeScan_resultsTables_simpleModel_11params_RECON2_newSSminDP10maxDP200_Dsyl_3RunsCombined_PLS20_Joint_retSims10000_gridPoints33_FINAL")

### Plotting parameters
# Define number of chromosomes
nChr <- 15
# Should you wish to plot a horizontal line representing the credible interval, we can set this here. Note, this DOES NOT define the significance states of the input table, rather this is only used for plotting purposes.
credible_interval <- 0.999
#credible_interval_1D <- round(credible_interval^0.5,3)
credible_interval_1D <- 0.995

# Cap maximum likelihood, in case some values are extremely high.
likelihood_max_cap_2D <- 6
likelihood_max_cap_1D <- 6
# Cap maximum and minumum asymmetry, in case some values are extremely high.
asymmetry_min_max_cap_value <- 3

# Parameter labels
param_labels <- c("migHL", "migLH")

# Input file
input_file <- "genomeScan_table.allChr1_15.txt"

### Plot

# Read in genome scan table
genomeScan_table <- read.delim(paste0(input_file), sep = " ")

# Check data classes
sapply(genomeScan_table, class)
# Rename headers (in case header special characters modified)
colnames(genomeScan_table) <- c("scaffold", "Chr", "Chr_position", "window_start", "window_end", "-log(1-p_migHL)", "-log(1-p_migLH)", "-log(1-p_joint_migHL_migLH)", paste0("significance_",credible_interval_1D,"_migHL"), paste0("significance_",credible_interval_1D,"_migLH"), paste0("significance_",credible_interval,"_joint_migHL_migLH"), "mode_migHL", "mode_migLH", "asymmetry")

# Cap max -log(1-p) to a reasonable value, e.g. 5
cap_max_prob <- function(x, likelihood_max_cap){
  if(x > likelihood_max_cap){
    x <- likelihood_max_cap - runif(1,0,1) # add jitter to distinguish between capped points
    return(x)
  }
  else {
    return(x)
  }
}
for (prob_column in c(6,7)) {
  genomeScan_table[,prob_column] <- sapply(genomeScan_table[,prob_column], cap_max_prob, likelihood_max_cap=likelihood_max_cap_1D)
}
for (prob_column in c(8)) {
  genomeScan_table[,prob_column] <- sapply(genomeScan_table[,prob_column], cap_max_prob, likelihood_max_cap=likelihood_max_cap_2D)
}

# Let's transform asymmetry values to emphasize extreme values, here via log(a/(1-a)), known as the odds ratio of the posterior
genomeScan_table[,ncol(genomeScan_table)+1] <- log(genomeScan_table$asymmetry / (1 - genomeScan_table$asymmetry))
colnames(genomeScan_table)[ncol(genomeScan_table)-1] <- "asymmetry_original"
colnames(genomeScan_table)[ncol(genomeScan_table)] <- "asymmetry"

# To color significant values by asymmetry
# Let's cap (new) asymmetry values above 3 to 3. Recall log(a/(1-a)) = 3 corresponds to a = 0.95.
cap_min_max_asymmetry <- function(x, asymmetry_min_max_cap){
  if(x > asymmetry_min_max_cap){
    x <- asymmetry_min_max_cap
    return(x)
  }
  else if(x < -asymmetry_min_max_cap){
    x <- -asymmetry_min_max_cap
    return(x)
  }
  else {
    return(x)
  }
}
genomeScan_table$asymmetry <- sapply(genomeScan_table$asymmetry, cap_min_max_asymmetry, asymmetry_min_max_cap=asymmetry_min_max_cap_value)

# To distinguish chromosomes, we plot chromosomes in alternating shades of gray. To do this, we define as numeric the "Chr" column, which results in numbers 1-15, and we return the modulus given 2, which returns either 0 or 1 (i.e. even or odd). Then we define this as factor for which to color.
genomeScan_table[,ncol(genomeScan_table)+1] <- as.data.frame(as.numeric(genomeScan_table$Chr))
genomeScan_table[,ncol(genomeScan_table)+1] <- as.factor(genomeScan_table[,ncol(genomeScan_table)]%%2)
# To color significant values by asymmetry
genomeScan_table[,ncol(genomeScan_table)+1] <- (genomeScan_table$significance_0.999_joint_migHL_migLH== 1)*genomeScan_table$asymmetry
genomeScan_table[,ncol(genomeScan_table)][genomeScan_table[,ncol(genomeScan_table)] == 0] <- NA
# Define number of colour breaks and colour palette
#Option 1: Equal sized bins across colour scale
# col_breaks <- 11
# #genomeScan_table[,ncol(genomeScan_table)+1] <- as.numeric(cut(genomeScan_table[,ncol(genomeScan_table)],breaks = col_breaks)) # breaks with min/max value = min/max of data
# #genomeScan_table[,ncol(genomeScan_table)+1] <- round(genomeScan_table[,ncol(genomeScan_table)]*10) + 1  # breaks with min/max value = 0/1 (fixed interval)
# genomeScan_table[,ncol(genomeScan_table)+1] <- round(genomeScan_table[,ncol(genomeScan_table)]*(5/3)) + 6  # breaks with min/max value = 0/1 (fixed interval)
# col_pal <- colorRampPalette(c('red3','white','dodgerblue3'))
# genomeScan_table[,ncol(genomeScan_table)+1] <- col_pal(col_breaks)[genomeScan_table[,ncol(genomeScan_table)]]
# Option 2: Equal sized bins across colour scale, apart from the central bin around asymmetry=0, which is larger
genomeScan_table[,ncol(genomeScan_table)+1] <- round(genomeScan_table[,ncol(genomeScan_table)]*(14/3)) + 15  # breaks with min/max value = 0/1 (fixed interval)
cols_asymmetry <- c("#CD0000", "#D01212", "#D42424", "#D73636", "#DB4848", "#DE5B5B", "#E26D6D", "#E67F7F", "#E99191", "#EDA3A3", "#F0B6B6", "#F4C8C8", "#F7DADA", "#FFFFFF", "#FFFFFF", "#FFFFFF", "#DEEBF7", "#CDE1F4", "#BDD7F0", "#ACCDED", "#9CC3E9", "#8BB9E6", "#7BAFE2", "#6AA5DE", "#5A9BDB", "#4991D7", "#3987D4", "#287DD0", "#1874CD")
genomeScan_table[,ncol(genomeScan_table)+1] <- cols_asymmetry[genomeScan_table[,ncol(genomeScan_table)]]
genomeScan_table[,ncol(genomeScan_table)+1] <- genomeScan_table[,ncol(genomeScan_table)]
genomeScan_table[,ncol(genomeScan_table)][!is.na(genomeScan_table[,ncol(genomeScan_table)])] <- "black"

# Let's get chromosome break points
chr_breaks <- c()
chr_names <- c()
for (i in seq(1:nChr)) {
  chr_break <- sort(which(genomeScan_table$Chr==paste0("Dsyl_",i,"_ml")))[1]
  chr_breaks <- c(chr_breaks,chr_break)
  chr_name <- paste0("Chr",i)
  chr_names <- c(chr_names,chr_name)
}
# Add final break (end) point
chr_break <- sort(which(genomeScan_table$Chr==paste0("Dsyl_",nChr,"_ml")), decreasing = TRUE)[1]
chr_breaks <- c(chr_breaks,chr_break)
# And now, get chromosome mid-points
chr_bins <- c()
for (j in seq(1:nChr)) {
  chr_bin <- (chr_breaks[j] + chr_breaks[j+1])/2
  chr_bins <- c(chr_bins,chr_bin)
}

# Manhatten plot based on joint 2D credible intervals
par(mfrow=c(1,1),mar=c(4,4,4,8))
# For alternate shades for dark and light grays
#plot(seq(1,nrow(genomeScan_table)), genomeScan_table$`-log(1-p_joint_migHL_migLH)`, xlab="Position", ylab="-log10(1-P)", main = "Joint 2D migHL-migLH significance - All Chromosomes", pch=19, cex = 0.25, col=c("gray0", "gray55")[genomeScan_table[,ncol(genomeScan_table)-3]])
plot(seq(1,nrow(genomeScan_table)), genomeScan_table$`-log(1-p_joint_migHL_migLH)`, xlab="Position", ylab="-log10(1-P)", main = "Joint 2D migHL-migLH significance - All Chromosomes", xaxt = "n", pch=19, cex = 0.25, col=c("gray0", "gray55")[genomeScan_table[,ncol(genomeScan_table)-4]])
axis(1, at=chr_bins, labels=chr_names)
# For distinguishing joint and marginal significance points
points(seq(1,nrow(genomeScan_table)), genomeScan_table$`-log(1-p_joint_migHL_migLH)`, pch=21, cex = 1.5, col=genomeScan_table[,ncol(genomeScan_table)], bg=genomeScan_table[,ncol(genomeScan_table)-1])
#points(seq(1,nrow(genomeScan_table)), genomeScan_table$`-log(1-p_joint_migHL_migLH)`, pch=19, cex = 0.25, col=c("red3", NA)[genomeScan_table$significance_0.999_joint_migHL_migLH])
#points(seq(1,nrow(genomeScan_table)), genomeScan_table$`-log(1-p_joint_migHL_migLH)`, pch=24, cex = 0.5, col=c(NA, NA)[genomeScan_table$significance_0.995_migHL], bg=c("red3", NA)[genomeScan_table$significance_0.995_migHL])
#points(seq(1,nrow(genomeScan_table)), genomeScan_table$`-log(1-p_joint_migHL_migLH)`, pch=25, cex = 0.5, col=c(NA, NA)[genomeScan_table$significance_0.995_migLH], bg=c("red3", NA)[genomeScan_table$significance_0.995_migLH])
# For unique colours per chromosome
#palette(c("darkgoldenrod2", "aquamarine2", "darkgray", "coral1", "dodgerblue3", "darkolivegreen2", "darkorchid4", "darkolivegreen4", "turquoise3", "firebrick2", "hotpink1", "wheat", "red4", "yellow2", "slateblue1"))
#plot(seq(1,nrow(genomeScan_table)), genomeScan_table$`-log(1-p_joint_migHL_migLH)`, xlab="Position", ylab="-log10(1-P)", main = "Joint 2D migHL-migLH significance - All Chromosomes", pch=19, cex = 0.3, col=genomeScan_table[,2])
abline(h = -log10(1-credible_interval), lty = 2, lwd = 1.5, col = "grey25")
#legend("bottomright", legend = c("marginal mHL", "marginal mLH", "marginal mHL + marginal mLH"), pt.bg=c("red3","red3", "red3"), col=c("red3","red3", "red3"), pch = c(24, 25, 11), bty = "o", pt.cex = 1, cex = 0.75, horiz = FALSE, inset = c(0.01, 0.025), bg = "grey95", text.col = "gray15", box.lwd = 0)
#ColorLegend("right", col=col_pal(col_breaks), labels=sprintf("%.1f",seq(0,1,0.1)), cntrlbl = TRUE, cex=0.5, inset = -0.015, width = round(nrow(genomeScan_table)/30), height = min(likelihood_max_cap_2D,max(genomeScan_table$`-log(1-p_joint_migHL_migLH)`)))
#ColorLegend("right", col=col_pal(col_breaks*5), labels=sprintf("%.1f",seq(-3,3,0.5)), adj = c(0,0.5), cntrlbl = TRUE, cex=0.5, inset = -0.04, width = round(nrow(genomeScan_table)*40), height = min(likelihood_max_cap_2D,max(genomeScan_table[,8]))*1.08)
ColorLegend("right", col=cols_asymmetry, labels=sprintf("%.1f",seq(-3,3,0.5)), cntrlbl = TRUE, cex=0.5, inset = -0.015, width = round(nrow(genomeScan_table)/30), height = min(likelihood_max_cap_2D,max(genomeScan_table$`-log(1-p_joint_migHL_migLH)`)))
#ColorLegend(x=54800000, y=8, col=col_pal(col_breaks), labels=sprintf("%.1f",seq(0,1,0.1)), cntrlbl = TRUE, cex=0.5, inset = -0.05, width = round(nrow(genomeScan_table)*40), height = 10)
#legend("topright", legend = c("asymmetry"), bty = "n", cex = 0.5, inset = 0)

# Manhatten plot based on 1D credible intervals
# par(mfrow=c(1,1),mar=c(4,4,4,4))
# for (p in 1:2) {
#   # We want to distinguish between windows that satisfy a single parameter CI signficance and those that satisfy the joint (2) parameter CI significance.
#   # We do this by adding a third value (3), which signifies windows that satisfy the 1D CI but not the joint 2D CI
#   signficance_combined_vector <- vector()
#   for (i in seq(1, nrow(genomeScan_table))) {
#     if ( (genomeScan_table[,p+8][i] == 1) & (genomeScan_table[,11][i] == 2)) {
#       signficance_combined_vector[i] <- 3
#     }
#     else {
#       signficance_combined_vector[i] <- genomeScan_table[,11][i]
#     }
#   }
#   # And plot
#   plot(seq(1,nrow(genomeScan_table)), genomeScan_table[,p+5], xlab="Position", ylab="-log10(1-P)", ylim = c(0, likelihood_max_cap_1D), main = paste0(param_labels[p], " (marginal) significance - All Chromosomes"), pch=19, cex = 0.25, col=c("gray0", "gray55")[genomeScan_table[,16]])
#   points(seq(1,nrow(genomeScan_table)), genomeScan_table[,p+5], pch=19, cex = 0.25, col=c("red3", NA,"darkorange")[signficance_combined_vector])
#   abline(h = -log10(1-credible_interval_1D), lty = 2, lwd = 1.5, col = "grey30")
#   #legend("topright", legend = c(paste0("Not significant(<SQRT(",credible_interval*100,")% CI)"), paste0("Significant for this marginal parameter (>SQRT(", credible_interval*100, ")% CI)"), paste0("Significant for (both) joint parameters (>", credible_interval*100,"% CI)")), col=c("grey35","darkorange", "red3"), pch = 19, bty = "n", pt.cex = 1.25, cex = 0.6, horiz = FALSE, inset = c(0.01, 0.025))
#   legend("topright", legend = c(paste0("Not significant(<",credible_interval_1D*100,"% CI)"), paste0("Significant for this marginal parameter (>", credible_interval_1D*100, "% CI)"), paste0("Significant for (both) joint parameters (>", credible_interval*100,"% CI)")), col=c("grey35","darkorange", "red3"), pch = 19, bty = "n", pt.cex = 1.25, cex = 0.6, horiz = FALSE, inset = c(0.01, 0.025))
# }
