library('emmeans')
library('car')
library('ggplot2')
library(glmmTMB)
library(ggeffects)
library(rgl)
library(patchwork)
library(DHARMa)

# This script tests fitnss effects of alternate CEN genotypes in parental plants (low elevation populations) growing in the transplant experiment (low elevation sites). 

lslg<-read.table('2.CEN_parentals_fitness_lslg.txt',header=T,sep='\t')
lslg2<-na.omit(lslg)
lslg2$Call[lslg2$Call == "Homozygous 1/1"] <- "LL"
lslg2$Call[lslg2$Call == "Homozygous 2/2"] <- "HH"


fit <- glmmTMB(seeds5years ~ Call * poly(ros5globYears,2) * poly(day5globYears,2), data=lslg2, family=nbinom2,ziformula =~.) 
# Explore models
fit_2a <- glmmTMB(seeds5years ~ Call * poly(ros5globYears,2) + Call * poly(day5globYears,2), data=lslg2, family=nbinom2,ziformula =~.)
fit_2c <- glmmTMB(seeds5years ~ Call * poly(day5globYears,2) + poly(ros5globYears,2), data=lslg2, family=nbinom2,ziformula =~.)
fit_2g <- glmmTMB(seeds5years ~ Call * poly(ros5globYears,2) * poly(day5globYears,2) + site + (1|site_plot), data=lslg2, family=nbinom2,ziformula =~.) 
#fit_2g <- glmmTMB(seeds5years ~ Call * poly(ros5globYears,2) * poly(day5globYears,2) + (1|site/site_plot), data=lslg2, family=nbinom2,ziformula =~.)  #  no convergence
fit_2h <- glmmTMB(seeds5years ~ Call * poly(ros5globYears,2) * poly(day5globYears,2) + site, data=lslg2, family=nbinom2,ziformula =~.)  
print(anova(fit_2,fit_2a))
print(anova(fit_2,fit_2c))
print(anova(fit_2,fit_2g))
print(anova(fit_2,fit_2h))
# fit is the best model

summary(fit)
plot(simulateResiduals(fit))
Anova(fit,type = 3, component="cond")
Anova(fit,type = 3, component="zi")

#write.table(Anova(fit,type = 3, component="cond"),'lslg_cond.txt',row.names = T,sep='\t',quote=F)
#write.table(Anova(fit,type = 3, component="zi"),'lslg_zi.txt',row.names = T,sep='\t',quote=F)

###### make plots ##########

dat<-ggpredict(fit, terms = c("day5globYears [all]","Call"), type="fe") 
plot2<-plot(dat) + theme(panel.grid = element_blank()) ; plot2 + ylim(0, 2000)
# Extract X-axis and colors
y_limits <- layer_scales(plot2)$x$range$range
layer_data <- layer_data(plot2)$colour
colors <- unique(layer_data)
# Make boxplot
boxplot <- ggplot(lslg, aes(x = Call, y = day5globYears, fill = Call)) +
  scale_fill_manual(values = colors) +
  scale_y_continuous(limits = y_limits) +
  geom_boxplot() +
  theme_bw() +
  coord_flip() +
  theme(panel.grid = element_blank()) +
  guides(fill = FALSE)
print(boxplot)
# Combine scatterplot and boxplot
combined_plot <- plot2 + boxplot + plot_layout(nrow = 2, heights = c(4, 1))
print(combined_plot)


