library(ggplot2)
library(lmerTest)
library(emmeans)

# This script plots and tests differences in flowering time of alternate CEN genotypes from parenatal genotypes growing in the reciprocal transplant experiment

df<-read.table('CEN_parentals_flowerintTime.csv',header=T,sep=',',na.string=c('NA','na','n.a.',''),fill=T)

# standardise FLOWERING TIME across 4 sites
df2<-df
df2$day16_globScaled<-scale(df2$day16)
df2$day17_globScaled<-scale(df2$day17)
df2$day18_globScaled<-scale(df2$day18)
df2$day19_globScaled<-scale(df2$day19)
df2$day20_globScaled<-scale(df2$day20)
df2$day5globYears<-rowMeans(df2[,c(19:23)],na.rm=T)

df2$Call <- factor(df2$Call, levels = c("Homozygous 1/1", "Heterozygous 1/2", "Homozygous 2/2"))  
df2$geno<-as.character(df2$Call)
df2$geno[df2$geno == "Homozygous 1/1"] <- "LL"
df2$geno[df2$geno == "Heterozygous 1/2"] <- "HL"
df2$geno[df2$geno == "Homozygous 2/2"] <- "HH"
df2$geno <- factor(df2$geno, levels = c("LL", "HL", "HH"))  
df2$ecotype<-as.character(df2$altitude)
df2$ecotype[df2$ecotype == "high"] <- "High ecotype"
df2$ecotype[df2$ecotype == "low"] <- "Low ecotype"
df2$ecotype <- factor(df2$ecotype, levels = c("High ecotype", "Low ecotype"))  
df2$trsite<-as.character(df2$site_altitude)
df2$trsite[df2$trsite == "high_site"] <- "High site"
df2$trsite[df2$trsite == "low_site"] <- "Low site"
df2$trsite <- factor(df2$trsite, levels = c("High site", "Low site")) 


# Supp. fig., inclduing both ecotypes
ggplot(df2, aes(x = geno, y = day5globYears, fill = geno)) +
  geom_boxplot() +
  geom_point(position = position_jitter(width = 0.2), color = "grey") +
  facet_grid(. ~ ecotype + trsite) +
  scale_fill_manual(values = c("red", "orange", "blue")) +
  labs(y = "Flowering time") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  guides(fill = FALSE)

# Main text, inclduing only low ecotype
lg<-subset(df2,altitude=='low')
ggplot(lg, aes(x = geno, y = day5globYears, fill = geno)) +
  geom_boxplot() +
  #geom_point(position = position_jitter(width = 0.2), color = "grey") +
  facet_grid(. ~ ecotype + trsite) +
  scale_fill_manual(values = c("red", "orange", "blue")) +
  labs(y = "Flowering time") +
  theme_bw() +
  theme(panel.grid = element_blank()) 

# Test differences and extract estimates
m3<-lmer(day5globYears ~ geno * site_altitude + altitude + (1|site/site_plot), data=df2, na.action=na.omit)
m3b<-lmer(day5globYears ~ geno + site_altitude + altitude + (1|site/site_plot), data=df2, na.action=na.omit)
anova(m3,m3b)
summary(m3)
emmeans(m3, pairwise~geno|site_altitude)


