library(lmerTest)
library(ggeffects)
library(ggplot2)
library(emmeans)
library(car)

# This script tests the effect of plant size on flowering  across sites and CEN genotypes

df<-read.table('4.CEN_F2_traits.txt',sep='\t',header=T)
df$CEN <- factor(df$CEN, levels = c("HL", "HH"))  # this is to have correct colors in plots automatically
colnames(df)[13:14]<-c('size','day')

# Test correlation between size and flowering time across sites and CEN genotypes
# Find best model:
fit0<-lmer(day ~ size * site * CEN * cluster + (1|site_plot), data=df, na.action=na.omit)
fit1<-lmer(day ~ size * site * CEN + cluster + (1|site_plot), data=df, na.action=na.omit)
fit2<-lmer(day ~ size * site + CEN + cluster + (1|site_plot), data=df, na.action=na.omit)
fit3<-lmer(day ~ size + site + CEN + cluster + (1|site_plot), data=df, na.action=na.omit)
fit4<-lmer(day ~ size + CEN + cluster + (1|site_plot), data=df, na.action=na.omit)  
fit5<-lmer(day ~ size + cluster + (1|site_plot), data=df, na.action=na.omit)
#test random effects
fit6<-lm(day ~ size + CEN + cluster, data=df, na.action=na.omit) # best model
print(anova(fit1,fit4))
print(anova(fit4,fit6)) 
Anova(fit6,type=3)
plot(fit6)
# plot
v<-ggpredict(fit6, terms = c("size [all]","CEN")) 
plot2<-plot(v) + geom_vline(xintercept=0, linetype="dashed"); plot2

emm_fit6 <- test(emmeans(fit6, "CEN"))
em_cont_fit6<-emmeans(fit6, specs = pairwise ~ CEN, type="linear", adjust="none")
em_cont_fit6
emm_fit6




