# This script calcualtes Dstatistics for the iterations and makes plots

# Requires:
# arg1: trio file
# arg2: collected iterations with adjusted sign (from previous script)
# arg3: plot name
# arg4: table name

# Rscript 7_calculate_Dstat.R   trios_BBAA.txt   collected_trios_adjustedSign_BBAA.txt   collected_trios_adjustedSign_BBAA_Dstat.pdf  quantile_table_BBAA.csv

args <- commandArgs(trailingOnly = TRUE)


df1<-read.table(args[1],header=T,sep='\t')
df1$trios<-paste(df1$P1,df1$P2,df1$P3,sep=',')

df2<-read.table(args[2],header=T,sep='\t')
df2$trios<-paste(df2$P1,df2$P2,df2$P3,sep=',')
df2$Dstatistics<-(df2$ABBA-df2$BABA)/(df2$ABBA+df2$BABA)

# Create an empty data frame to store the quantiles
quantile_table <- data.frame(
  trio = character(),
  q025 = numeric(),
  q975 = numeric(),
  meanD = numeric(),
  stringsAsFactors = FALSE )

# Loop trios
pdf(args[3], width = 12, height = 15)  
par(mfrow = c(5, 4))
for (t in unique(df1$trios)) {
  data <- subset(df2, trios == t)
  q25 <- quantile(data$Dstatistic, c(0.025, 0.975))
  meanD <- mean(data$Dstatistic)
  
  # Save quantiles
  quantile_table <- rbind(quantile_table, data.frame(
    trio = t,
    q025 = q25[1],
    q975 = q25[2],
    meanD = meanD
  ))
  
  # Title color logic
  title_color <- "black"
  if (q25[1] > 0 | q25[2] < 0) {
    title_color <- "red"
  }
  
  hist(data$Dstatistic, main = t, xlim = c(-1, 1), breaks = 50, 
       col = "darkgrey", border = "white", col.main = title_color)
  abline(v = meanD, col = "red", lwd = 2, lty = 1)
  abline(v = 0, col = "black", lwd = 1)
}

dev.off()

print("This is the quantile table:")
print(quantile_table)
write.csv(quantile_table, args[4], row.names = FALSE)



