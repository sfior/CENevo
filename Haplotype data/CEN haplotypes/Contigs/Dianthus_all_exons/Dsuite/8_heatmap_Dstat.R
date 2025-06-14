# Make your file 'collected_trios_adjustedSign_4Rplot.txt' with mean values from the quantile table
# Run Rscript 8_heatmap_Dstat.R 

library(ggplot2)
library(reshape2)
library(dplyr)

df <- read.table("collected_trios_adjustedSign_4Rplot.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Add reverse pairs to make the matrix symmetric
df_rev <- df %>% 
  rename(Taxon1 = Taxon2, Taxon2 = Taxon1)

df_full <- bind_rows(df, df_rev)

taxa_order <- c("Group_6", "Group_5", "Group_4", "Ds", "Group_3", "Group_2")

# Add diagonal entries (optional NA)
diag_df <- data.frame(Taxon1 = taxa_order, Taxon2 = taxa_order, D_BBAA = NA)
df_full <- bind_rows(df_full, diag_df)

# Create wide matrix
heatmap_matrix <- dcast(df_full, Taxon1 ~ Taxon2, value.var = "D_BBAA")
rownames(heatmap_matrix) <- heatmap_matrix$Taxon1
heatmap_matrix <- heatmap_matrix[, -1]

heatmap_matrix <- heatmap_matrix[taxa_order, taxa_order]

heatmap_long <- melt(as.matrix(heatmap_matrix), varnames = c("Taxon1", "Taxon2"), value.name = "D_BBAA")

heatmap_long$Taxon1 <- factor(heatmap_long$Taxon1, levels = taxa_order)
heatmap_long$Taxon2 <- factor(heatmap_long$Taxon2, levels = taxa_order)


pdf("Rplot_allDs_BBAA_grey.pdf", width = 6, height = 6)
ggplot(heatmap_long, aes(Taxon1, Taxon2, fill = D_BBAA)) +
  geom_tile(color = "grey50") +  
  scale_fill_gradient(
    low = "white", high = "black",  
    na.value = "burlywood1",        
    limits = c(min(heatmap_long$D_BBAA, na.rm = TRUE), max(heatmap_long$D_BBAA, na.rm = TRUE))
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_line(color = "grey80", size = 0.2)
  ) +
  coord_fixed() +
  labs(title = "D_BBAA Heatmap (White to Black)", fill = "D_BBAA")
dev.off()

