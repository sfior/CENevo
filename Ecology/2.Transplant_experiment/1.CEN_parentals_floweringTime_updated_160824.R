library(ggplot2)
library(lmerTest)
library(emmeans)
library(partR2)
library(dplyr)
# This script plots and tests differences in flowering time of alternate CEN genotypes from parental genotypes growing in the reciprocal transplant experiment

df<-read.table('1.CEN_parentals_flowerintTime.csv',header=T,sep=',',na.string=c('NA','na','n.a.',''),fill=T)

#Raw mean day for plotting
df$raw_mean_day <- rowMeans(df[, c('day16', 'day17', 'day18', 'day19', 'day20')], na.rm = TRUE)

#Standardise flowering time within elevational environment 
df2 <- df %>%
  group_by(site_altitude) %>%
  mutate(day16_Scaled = scale(day16),
         day17_Scaled = scale(day17),
         day18_Scaled = scale(day18),
         day19_Scaled = scale(day19),
         day20_Scaled = scale(day20)) %>%
  ungroup()

df2$day5ScaledYears <- rowMeans(df2[, c('day16_Scaled', 'day17_Scaled', 'day18_Scaled', 'day19_Scaled', 'day20_Scaled')], na.rm = TRUE)

#
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
ggplot(df2, aes(x = geno, y = raw_mean_day, fill = geno)) +
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
ggplot(lg, aes(x = geno, y = raw_mean_day, fill = geno)) +
  geom_boxplot() +
  #geom_point(position = position_jitter(width = 0.2), color = "grey") +
  facet_grid(. ~ ecotype + trsite) +
  scale_fill_manual(values = c("red", "orange", "blue")) +
  labs(y = "Flowering time") +
  theme_bw() +
  theme(panel.grid = element_blank()) 

#Test differences and extract estimates
df2 <- df2[!is.na(df2$day5ScaledYears), ]
m3<-lmer(day5ScaledYears ~ geno * site_altitude + altitude + (1|site/site_plot), data=df2, na.action=na.omit)
m3b<-lmer(day5ScaledYears ~ geno + site_altitude + altitude + (1|site/site_plot),data=df2, na.action=na.omit)
anova(m3,m3b)
summary(m3)
emmeans(m3, pairwise~geno|site_altitude)

#Backtransform estimates and contrasts for reporting of results
mean_sd_by_altitude <- df %>%
  group_by(site_altitude) %>%
  summarise(across(c(day16, day17, day18, day19, day20), 
                   list(mean = ~mean(.x, na.rm = TRUE), 
                        sd = ~sd(.x, na.rm = TRUE)), 
                   .names = "{col}_{fn}"))

emm_df_contrasts <- as.data.frame(emm_results$contrasts)

#Contrasts
emm_df_contrasts_high <- subset(emm_df_contrasts, site_altitude=="high_site")
contrasts_df_high <- emm_df_contrasts_high %>%
  dplyr::mutate(backtransformed_estimate = estimate * overall_sd_high)
contrasts_df_high <- contrasts_df_high %>%
  dplyr::mutate(backtransformed_SE = SE * overall_sd_high)

emm_df_contrasts_low <- subset(emm_df_contrasts, site_altitude=="low_site")
contrasts_df_low <- emm_df_contrasts_low %>%
  dplyr::mutate(backtransformed_estimate = estimate * overall_sd_low)
contrasts_df_low <- contrasts_df_low %>%
  dplyr::mutate(backtransformed_SE = SE * overall_sd_low)

emm_df_contrasts_backtransformed_day <- rbind(contrasts_df_high, contrasts_df_low)


#Proportion of variance explained full data
m3<-lmer(day5ScaledYears ~ geno * site_altitude+altitude + (1|site/site_plot), data=df2, na.action=na.omit)
R2_GPc_part1 <- partR2(m3, partvars = c("geno:site_altitude"), data = df2, 
                       nboot = 10000)

m3_1<-lmer(day5ScaledYears ~ geno +site_altitude+altitude + (1|site/site_plot), data=df2, na.action=na.omit)

R2_GPc_part2 <- partR2(m3_1, partvars = c("geno", "site_altitude", "altitude"), data = df2, 
                       nboot = 10000)
effect_of_geno_full_data<-mergeR2(R2_GPc_part1, R2_GPc_part2)


#Proportion of variance explain low ecotype 
low_eco <- subset(df2, altitude=="low")
m3<-lmer(day5ScaledYears ~ geno * site_altitude + (1|site/site_plot), data=low_eco, na.action=na.omit)
#r.squaredGLMM(m3)
summary(m3)
R2_GPc_part1 <- partR2(m3, partvars = c("geno:site_altitude"), data = low_eco, 
                       nboot = 10000)
m3_1<-lmer(day5ScaledYears ~ geno +site_altitude + (1|site/site_plot), data=low_eco, na.action=na.omit)
summary(m3)
R2_GPc_part2 <- partR2(m3_1, partvars = c("geno", "site_altitude"), data = low_eco, 
                       nboot = 10000)
effect_of_cen_low_eco_only <- mergeR2(R2_GPc_part1, R2_GPc_part2) 
