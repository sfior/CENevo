library(lmerTest)
library(ggeffects)
library(ggplot2)
library(visreg)
library(emmeans)
library(rgl)
library(car)
library(mgcv)

# This script performs analyses of correlational selection on alternate CEN genotypes growing in the transplant experiment

df<-read.table('5.CEN_F2_fitness.txt',sep='\t',header=T)
df$CEN <- factor(df$CEN, levels = c("HL", "HH"))  # this is to have correct colors in plots automatically
colnames(df)[13:14]<-c('size','day')

# correlational selection as in Svensson 2021
# Findeln reduced model
dfF<-subset(df, site == 'Findeln')
mF<-lmer(aster_fitness ~ CEN + size + day + 
               I(0.5 * size^2) + I(0.5 * day^2) + size:day + 
               CEN:size + CEN:day + CEN:size:day +
               I(0.5 * size^2):CEN + (1|cluster), data=dfF, na.action=na.omit)  
Anova(mF,type = 3)
plot(mF)
qqnorm(resid(mF))
qqline(resid(mF), col = "red")
hist(resid(mF), breaks = 20)

# Zeneggen reduced model
dfZ<-subset(df, site == 'Zeneggen')
mZ<-lmer(aster_fitness ~ CEN + size + day + 
               I(0.5 * day^2) + size:day + 
               CEN:size + CEN:day + CEN:size:day +
               I(0.5 * day^2):CEN + (1|cluster), data=dfZ, na.action=na.omit)  
Anova(mZ,type = 3)
plot(mZ)
qqnorm(resid(mZ))
qqline(resid(mZ), col = "red")
hist(resid(mZ), breaks = 20)

# Predict fitenss surfaces for high (Findeln) and low (Zeneggen) site
#Findeln:
v1<-seq(min(dfF$day,na.rm=T), max(dfF$day,na.rm=T),0.05)
v2<-seq(min(dfF$size,na.rm=T), max(dfF$size,na.rm=T),0.05)
v_Findeln<-ggpredict(mF, terms = c("day [v1]", "size [v2]", "CEN"))

v_Findeln_HH<-subset(v_Findeln, facet=='HH')
v_Findeln_HL<-subset(v_Findeln, facet=='HL')
x1 = as.numeric(as.character(v_Findeln_HL$x))
y1 = as.numeric(as.character(v_Findeln_HL$group))
z1 = as.numeric(as.character(v_Findeln_HL$predicted))
x2 = as.numeric(as.character(v_Findeln_HH$x))
y2 = as.numeric(as.character(v_Findeln_HH$group))
z2 = as.numeric(as.character(v_Findeln_HH$predicted))

# For analyses separate by site:
#Zeneggen:
v1<-seq(min(dfZ$day,na.rm=T), max(dfZ$day,na.rm=T),0.075)
v2<-seq(min(dfZ$size,na.rm=T), max(dfZ$size,na.rm=T),0.075)
v_Zeneggen<-ggpredict(mZ, terms = c("day [v1]", "size [v2]", "CEN"))

v_Zeneggen_HH<-subset(v_Zeneggen, facet=='HH')
v_Zeneggen_HL<-subset(v_Zeneggen, facet=='HL')
x3 = as.numeric(as.character(v_Zeneggen_HL$x))
y3 = as.numeric(as.character(v_Zeneggen_HL$group))
z3 = as.numeric(as.character(v_Zeneggen_HL$predicted))
x4 = as.numeric(as.character(v_Zeneggen_HH$x))
y4 = as.numeric(as.character(v_Zeneggen_HH$group))
z4 = as.numeric(as.character(v_Zeneggen_HH$predicted))

# Make 3D plots
mfrow3d(1, 2, sharedMouse = TRUE)
par3d(userMatrix = NULL)
plot3d(x1, y1, z1, type = "p", col = "red", size = 2, xlab = 'Flowering time', ylab='Size',zlab='Fitness')
plot3d(x2, y2, z2, type = "p", col = "blue", size = 2, add=T)

bottom_plane1 <- rbind(c(min(x1), min(y1), min(z1)), 
                       c(min(x1), median(y1), min(z1)),
                       c(max(x1), median(y1), min(z1)),
                       c(max(x1), min(y1), min(z1)))
bottom_plane2 <- rbind(c(min(x1), min(y1), min(z1)), 
                       c(min(x1), median(y1), min(z1)),
                       c(median(x1), median(y1), min(z1)),
                       c(median(x1), min(y1), min(z1)))
top_plane2 <- rbind(c(min(x1), min(y1), max(z2)),   
                    c(min(x1), median(y1), max(z2)),
                    c(median(x1), median(y1), max(z2)),
                    c(median(x1), min(y1), max(z2)))
polygon3d(top_plane2, col = 'lightgray', alpha=0.2)
myview <- readRDS("myview.rds")
par3d(userMatrix = myview)

next3d()
plot3d(x3, y3, z3, type = "p", col = "red", size = 2, xlab = 'Flowering time', ylab='Size',zlab='Fitness')
plot3d(x4, y4, z4, type = "p", col = "blue", size = 2, add=T)

bottom_plane3 <- rbind(c(min(x3), max(y3), min(z3)), 
                       c(max(x3), max(y3), min(z3)),
                       c(max(x3), median(y3), min(z3)),
                       c(min(x3), median(y3), min(z3)))
bottom_plane4 <- rbind(c(min(x3), max(y3), min(z3)), 
                       c(median(x3), max(y3), min(z3)),
                       c(median(x3), median(y3), min(z3)),
                       c(min(x3), median(y3), min(z3)))
top_plane4 <- rbind(c(min(x3), max(y3), max(z3)), 
                    c(median(x3), max(y3), max(z3)),
                    c(median(x3), median(y3), max(z3)),
                    c(min(x3), median(y3), max(z3)))
polygon3d(top_plane4, col = 'lightgray', alpha=0.2)
par3d(userMatrix = myview)


# Compare fitness at center of the polygons. Non-overlapping CI indicate signicant difference
# Findeln
x1<-min(dfF$size,na.rm=T)+( (max(dfF$size,na.rm=T) - min(dfF$size,na.rm=T)) /4)
y1<-min(dfF$day,na.rm=T)+( (max(dfF$day,na.rm=T) - min(dfF$day,na.rm=T)) /4)
vF<-ggpredict(mF, terms = c("day [y1]", "size [x1]", "CEN")) 
vF
# Zeneggen
y2<-min(dfZ$day,na.rm=T)+( (max(dfZ$day,na.rm=T) - min(dfZ$day,na.rm=T)) /4)   
x2<-max(dfZ$size,na.rm=T)-( (max(dfZ$size,na.rm=T) - min(dfZ$size,na.rm=T)) /4)
vZ<-ggpredict(mZ, terms = c("day [y2]", "size [x2]", "CEN")) 
vZ

### Fit with GAM
# Findeln
gamF <- gam(aster_fitness ~ 
                   te(size, day, by = CEN) +  # Tensor
                   s(size, by = CEN) +  # smooth size
                   s(day, by = CEN) +  # smooth day
                   I(size^2) +  # quadratic size
                   I(day^2) +  # quadratic day
                   CEN,
                 random = list(cluster = ~1),  # Random effect for cluster
                 data = dfF
)
summary(gamF)
gam.check(gamF)  
qq.gam(gamF, cex=0.5, pch=16)

# Zeneggen
gamZ <- gam(aster_fitness ~ 
              te(size, day, by = CEN) +  # Tensor
              s(size, by = CEN) +  # smooth size
              s(day, by = CEN) +  # smooth day
              I(size^2) +  # quadratic size
              I(day^2) +  # quadratic day
              CEN,
            random = list(cluster = ~1),  # Random effect for cluster
            data = dfZ
)
summary(gamZ)
gam.check(gamZ)  
qq.gam(gamZ, cex=0.5, pch=16)

# 3D plots
#dfF$CEN <- as.factor(dfF$CEN)
cen_levels <- levels(dfF$CEN)

mfrow3d(1, 2, sharedMouse = TRUE)
par3d(userMatrix = NULL)
# Findeln
x_vals <- seq(min(dfF$day, na.rm = TRUE), 
              max(dfF$day, na.rm = TRUE), length = 100)  # Switch: 'day' on x-axis
y_vals <- seq(min(dfF$size, na.rm = TRUE), 
              max(dfF$size, na.rm = TRUE), length = 100)  # Switch: 'size' on y-axis
grid <- expand.grid(day = x_vals, size = y_vals, CEN = cen_levels)  # Match switched variables
grid$pred <- predict(gamF, newdata = grid, type = "response")
colors <- c("HH" = "blue", "HL" = "red")
for (cen in cen_levels) {
  subset_grid <- subset(grid, CEN == cen)
  z_matrix <- matrix(subset_grid$pred, nrow = length(x_vals), ncol = length(y_vals))
  persp3d(x = x_vals, y = y_vals, z = z_matrix,
          col = colors[cen], alpha = 0.7, add = ifelse(cen == cen_levels[1], FALSE, TRUE),
          xlab = "Flowering time",
          ylab = "Size",
          zlab = "Fitness")
}
par3d(userMatrix = myview)
next3d()
# Zeneggen
x_vals <- seq(min(dfZ$day, na.rm = TRUE), 
              max(dfZ$day, na.rm = TRUE), length = 100)  # Switch: 'day' on x-axis
y_vals <- seq(min(dfZ$size, na.rm = TRUE), 
              max(dfZ$size, na.rm = TRUE), length = 100)  # Switch: 'size' on y-axis
grid <- expand.grid(day = x_vals, size = y_vals, CEN = cen_levels)  # Match switched variables
grid$pred <- predict(gamZ, newdata = grid, type = "response")
colors <- c("HH" = "blue", "HL" = "red")
for (cen in cen_levels) {
  subset_grid <- subset(grid, CEN == cen)
  z_matrix <- matrix(subset_grid$pred, nrow = length(x_vals), ncol = length(y_vals))
  persp3d(x = x_vals, y = y_vals, z = z_matrix,
          col = colors[cen], alpha = 0.7, add = ifelse(cen == cen_levels[1], FALSE, TRUE),
          xlab = "Flowering time",
          ylab = "Size",
          zlab = "Fitness")
}
par3d(userMatrix = myview)



