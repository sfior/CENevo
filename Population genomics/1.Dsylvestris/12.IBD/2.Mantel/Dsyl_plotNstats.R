library(raster)
library(RColorBrewer)
library(ade4)

# This script runs a Mantel test to test IDB (IBD vs Resistance) and plots the Circuitscape current map


dmatrix<-read.table(file='Output_resistances.txt',sep=" ",header=T,row.names=1)
gmatrix<-read.table(file='genetic_distances_Fst.txt',sep="\t",header=T)
gmatrix2 <- as.data.frame(lapply(gmatrix, function(x) x / (1 - x)))
gmatrix <- gmatrix2
# Perform Mantel test
m<-mantel.rtest(as.dist(dmatrix), as.dist(gmatrix), nrepet = 9999)

#Extract the lower triangular
d <- dmatrix[lower.tri(dmatrix)]
g <- gmatrix[lower.tri(gmatrix)]

# Combine into a data frame for plotting
xy <- data.frame(d = d, g = g, pair = 1:length(g))
# Plot the scatterplot
plot(xy$d, xy$g, xlab = 'distance (m)', ylab = 'Fst/(1-Fst)', pch = 16, cex = 2, main = 'IBD')
# Fit and plot the regression line
fit <- lm(g ~ d, data = xy)
abline(fit)
# Spearman correlation test
spearman <- cor.test(xy$d, xy$g, method = "spearman", alternative = "g", exact = FALSE)
# Add a legend with stats
legend("topleft", bty = "n", legend = paste(
  'y =', format(coef(fit)[[2]], digits = 4), 'x', '+', format(coef(fit)[[1]], digits = 4), "\n",
  "R2 = ", format(summary(fit)$adj.r.squared, digits = 4), "\n",
  "Spearman rho = ", format(spearman$estimate, digits = 4), "\n",
  "pvalue = ", format(spearman$p.value, digits = 4), "\n",
  "\n", "pvalue mantel test = ", m$pvalue, sep = ""
))
# Label points
text(xy$d, xy$g + 0.001, xy$pair)


# Plot Circuitscape current
e_max <- extent(7, 8, 45.5, 46.5)
ELEV.GMTED <- raster("Output_cum_curmap.asc")
ELEV.GMTED.Wallis <- crop(ELEV.GMTED, e_max)
plot(ELEV.GMTED.Wallis)
# check resolution of raster (expressed in lat and long degree)
res(ELEV.GMTED)






