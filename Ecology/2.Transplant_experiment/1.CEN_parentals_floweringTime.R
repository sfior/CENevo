rm(list = ls())
gc()
library(ggplot2)
library(lmerTest)
library(emmeans)
library(partR2)
library(dplyr)
library(MuMIn)
library(car)
# This script plots and tests differences in flowering time of alternate CEN genotypes from parenatal genotypes growing in the reciprocal transplant experiment
df<-read.table('1.CEN_parentals_flowerintTime.csv',header=T,sep=',',na.string=c('NA','na','n.a.',''),fill=T)

df$raw_mean_day <- rowMeans(df[, c('day16', 'day17', 'day18', 'day19', 'day20')], na.rm = TRUE)
# standardise flowering time with elevational environment 
df2 <- df %>%
  group_by(site_altitude) %>%
  mutate(day16_Scaled = scale(day16),
         day17_Scaled = scale(day17),
         day18_Scaled = scale(day18),
         day19_Scaled = scale(day19),
         day20_Scaled = scale(day20)) %>%
  ungroup()

df2$day5ScaledYears <- rowMeans(df2[, c('day16_Scaled', 'day17_Scaled', 'day18_Scaled', 'day19_Scaled', 'day20_Scaled')], na.rm = TRUE)


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

# Test differences and extract estimates
df2 <- df2[!is.na(df2$day5ScaledYears), ]
df2$site_altitude<-as.factor(df2$site_altitude)
df2$site<-as.factor(df2$site)
df2$site_plot<-as.factor(df2$site_plot)
df2$altitude<-as.factor(df2$altitude)

m1<-lmer(day5ScaledYears ~ geno * site_altitude * altitude + (1|site/site_plot), data=df2, na.action=na.omit)
Anova(m1)

m2<-lmer(day5ScaledYears ~ geno * site_altitude + altitude + (1|site/site_plot), data=df2, na.action=na.omit)
m2b<-lmer(day5ScaledYears ~ geno + site_altitude + altitude + (1|site/site_plot),data=df2, na.action=na.omit)
anova(m2,m2b)
summary(m2)
emmeans(m2, pairwise~geno|site_altitude)

#backtransform estimates and contrasts for reporting of results
mean_sd_by_altitude <- df %>%
  group_by(site_altitude) %>%
  summarise(across(c(day16, day17, day18, day19, day20), 
                   list(mean = ~mean(.x, na.rm = TRUE), 
                        sd = ~sd(.x, na.rm = TRUE)), 
                   .names = "{col}_{fn}"))

emm_results <- emmeans(m2, pairwise ~ geno | site_altitude)
emm_df_means <- as.data.frame(emm_results$emmeans)
emm_df_contrasts <- as.data.frame(emm_results$contrasts)
#write.csv(emm_df_contrasts, file="emm_df_contrasts.csv")
#Means
high_sites_mean_sd <- mean_sd_by_altitude[1, ]
overall_mean_high <- mean(unlist(high_sites_mean_sd[1, c("day16_mean", "day17_mean", "day18_mean", "day19_mean", "day20_mean")]))
overall_sd_high <- mean(unlist(high_sites_mean_sd[1, c("day16_sd", "day17_sd", "day18_sd", "day19_sd", "day20_sd")]))
emm_df_means_high<-subset(emm_df_means,site_altitude=="high_site")
emm_df_means_high$backtransformed_days <- emm_df_means_high$emmean * overall_sd_high + overall_mean_high
emm_df_means_high$backtransformed_lower_CL <- emm_df_means_high$lower.CL * overall_sd_high + overall_mean_high
emm_df_means_high$backtransformed_upper_CL <- emm_df_means_high$upper.CL * overall_sd_high + overall_mean_high

low_sites_mean_sd <- mean_sd_by_altitude[2, ]
overall_mean_low <- mean(unlist(low_sites_mean_sd[1, c("day16_mean", "day17_mean", "day18_mean", "day19_mean", "day20_mean")]))
overall_sd_low <- mean(unlist(low_sites_mean_sd[1, c("day16_sd", "day17_sd", "day18_sd", "day19_sd", "day20_sd")]))
emm_df_means_low<-subset(emm_df_means,site_altitude=="low_site")
emm_df_means_low$backtransformed_days <- emm_df_means_low$emmean * overall_sd_low + overall_mean_low
emm_df_means_low$backtransformed_lower_CL <- emm_df_means_low$lower.CL * overall_sd_low + overall_mean_low
emm_df_means_low$backtransformed_upper_CL <- emm_df_means_low$upper.CL * overall_sd_low + overall_mean_low

emm_df_means_backtransformed_day <- rbind(emm_df_means_high, emm_df_means_low)
emm_df_means_backtransformed_day

#Contrasts
emm_df_contrasts_high <- subset(emm_df_contrasts, site_altitude=="high_site")
contrasts_df_high <- emm_df_contrasts_high %>%
  dplyr::mutate(backtransformed_estimate = estimate * overall_sd_high)
contrasts_df_high <- contrasts_df_high %>%
  dplyr::mutate(backtransformed_SE = SE * overall_sd_high)

#Contrasts
emm_df_contrasts_low <- subset(emm_df_contrasts, site_altitude=="low_site")
contrasts_df_low <- emm_df_contrasts_low %>%
  dplyr::mutate(backtransformed_estimate = estimate * overall_sd_low)
contrasts_df_low <- contrasts_df_low %>%
  dplyr::mutate(backtransformed_SE = SE * overall_sd_low)

emm_df_contrasts_backtransformed_day <- rbind(contrasts_df_high, contrasts_df_low)

emm_df_contrasts_backtransformed_day
#write.csv(emm_df_contrasts_backtransformed_day, file="emm_df_contrasts_backtransformed_day.csv")


#Proportion of variance explain low ecotype 
low_eco <- subset(df2, altitude=="low")
m3<-lmer(day5ScaledYears ~ geno * site_altitude + (1|site/site_plot), data=low_eco, na.action=na.omit)
r.squaredGLMM(m3)
summary(m3)
R2_GPc_part1 <- partR2(m3, partvars = c("geno:site_altitude"), data = low_eco, 
                       nboot = 10000)

m3_1<-lmer(day5ScaledYears ~ geno +site_altitude + (1|site/site_plot), data=low_eco, na.action=na.omit)
summary(m3)
R2_GPc_part2 <- partR2(m3_1, partvars = c("geno", "site_altitude"), data = low_eco, 
                       nboot = 10000)
effect_of_cen_low_eco_only <- mergeR2(R2_GPc_part1, R2_GPc_part2) 

effect_of_cen_low_eco_only1<- as.data.frame(effect_of_cen_low_eco_only$R2)
#write.csv(effect_of_cen_low_eco_only1, file="effect_of_cen_low_data.csv")

