library(ggbiplot)
library(lmerTest)
library(emmeans)

##### PCA
# Dsylvestris
eschikon.mat<-read.table('01_Dsyl_phenoTraits.txt',sep='\t',header=T,na.string='NA')
traits<-c(3:17)
eschikon.pops<-eschikon.mat[,1]
eschikon.pops <- factor(eschikon.pops, levels = c("Saxon", "Saviese", "Varen","Tsanfleuron", "Challer", "Mauvoisin"))
eschikon.traits<-eschikon.mat[,traits]
eschikon.pca <- prcomp(eschikon.traits,center=TRUE,scale.=TRUE)

g <- ggbiplot(eschikon.pca,obs.scale=1,var.scale=1,groups=eschikon.pops,ellipse=F,circle=F,var.axes=F)
g <- g + scale_color_manual(name = 'Populations',values=c("red","orange","yellow3","blue","purple","turquoise"))
g <- g + scale_shape_manual(name="Populations",values=c(16,16,16,17,17,17))
g <- g + geom_point(aes(colour=eschikon.pops,shape=eschikon.pops),size=3)
g <- g + theme(legend.direction='horizontal',legend.position='top')+theme_bw()
g <- g + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
g <- g + theme(aspect.ratio=1)
print(g)

g <- ggbiplot(eschikon.pca, obs.scale = 1, var.scale = 1, ellipse = F, circle = T, labels =  eschikon.mat[,1], labels.size = 0) 
g <- g + scale_color_discrete(name = '') + theme(legend.direction = 'horizontal', legend.position = 'top') + theme_bw() 
g <- g + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
g <- g + theme(aspect.ratio=1)
print(g)


#Dcarthusianorum
eschikon.mat<-read.table('01_Dcar_phenoTraits.txt',sep='\t',header=T,na.string='NA')
traits<-c(3:17)
eschikon.pops<-eschikon.mat[,1]
eschikon.pops <- factor(eschikon.pops, levels = c("Niedergampel", "Unterstalden", "Grengiols","Faldumalp", "Gibidumsee", "Simplonpass"))
eschikon.traits<-eschikon.mat[,traits]
eschikon.pca <- prcomp(eschikon.traits,center=TRUE,scale.=TRUE)

g <- ggbiplot(eschikon.pca,obs.scale=1,var.scale=1,groups=eschikon.pops,ellipse=F,circle=F,var.axes=F)
g <- g + scale_color_manual(name = 'Populations',values=c("red","orange","yellow3","blue","purple","turquoise"))
g <- g + scale_shape_manual(name="Populations",values=c(16,16,16,17,17,17))
g <- g + geom_point(aes(colour=eschikon.pops,shape=eschikon.pops),size=3)
g <- g + theme(legend.direction='horizontal',legend.position='top')+theme_bw()
g <- g + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
g <- g + theme(aspect.ratio=1)
print(g)

g <- ggbiplot(eschikon.pca, obs.scale = 1, var.scale = 1, ellipse = F, circle = T, labels =  eschikon.mat[,1], labels.size = 0) 
g <- g + scale_color_discrete(name = '') + theme(legend.direction = 'horizontal', legend.position = 'top') + theme_bw() 
g <- g + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
g <- g + theme(aspect.ratio=1)
print(g)


# combined Dsylvestris + Dcarthusian
data1<-read.table('01_Dsyl_phenoTraits.txt',sep='\t',header=T,na.string='NA')
data2<-read.table('01_Dcar_phenoTraits.txt',sep='\t',header=T,na.string='NA')
# make Z scores
data1Z<-data1
data1Z[3:17]<-scale(data1Z[3:17])
data2Z<-data2
data2Z[3:17]<-scale(data2Z[3:17])
mergedZ<-rbind(data1Z,data2Z)

mergedZ.pca <- prcomp(mergedZ[3:17],center=TRUE,scale.=TRUE)

g <- ggbiplot(mergedZ.pca,obs.scale=1,var.scale=1,groups=mergedZ$Ecotype,ellipse=T,circle=F,var.axes=F,alpha=0)
g <- g + scale_color_manual(name = 'Ecotype',values=c(5,6,"blue","red"))
g <- g + scale_shape_manual(name="Ecotype",values=c(17,16,17,16))
g <- g + geom_point(aes(colour=mergedZ$Ecotype,shape=mergedZ$Ecotype)) 
g <- g + theme(legend.direction='horizontal',legend.position='top')+theme_bw() 
g <- g + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
g <- g + theme(aspect.ratio=1)
print(g)

######  Boxplots and statistics flowering time
### Dsylvestris
eschikon.mat<-read.table('01_Dsyl_phenoTraits.txt',sep='\t',header=T,na.string='NA')
eschikon.mat$Population <- factor(eschikon.mat$Population, levels = c("Saxon","Varen","Saviese", "Tsanfleuron", "Challer", "Mauvoisin"))
# log transform
eschikon.mat[,3]<-log(as.numeric(eschikon.mat[,3]))
# remove 2SD oultiers
tmp1<-cbind(eschikon.mat[,c(1,3)])
y<-c()
pops<-c("Saxon","Varen","Saviese", "Tsanfleuron", "Challer", "Mauvoisin") # must be in the same order as in table
for (p in pops){
  tmp2<-subset(tmp1,Population==p)
  x<-tmp2[,2]
  x[x > (mean(x,na.rm=T)+(2*sd(x,na.rm=T)))] <- NA
  x[x < (mean(x,na.rm=T)-(2*sd(x,na.rm=T)))] <- NA
  y<-c(y,x)
}
mat <- data.frame(
  Population = eschikon.mat$Population,
  Flowering_date.d. = y
)
mat$Population <- factor(mat$Population, levels = c("Saxon","Varen","Saviese", "Tsanfleuron", "Challer", "Mauvoisin"))
# make plots
boxplot(Flowering_date.d. ~ Population, data = mat, col = c('red', 'red', 'red', 'blue', 'blue', 'blue'),las=2)
#Test difference in flowering time between high- and low-elevation
mat<-data.frame(cbind(mat,eschikon.mat$Ecotype,eschikon.mat$Motherline))
colnames(mat)[3:4]<-c('Ecotype','Motherline')
mat$Flowering_date.d.<-as.numeric(as.character(mat$Flowering_date.d.))
m1<-lmer(Flowering_date.d. ~ Ecotype + (1|Population/Motherline), data=mat,  na.action=na.exclude)
m2 <- lmer(Flowering_date.d. ~ 1 + (1|Population/Motherline), data=mat,  na.action=na.exclude)
anova(m1,m2)
lsmeans(m1, 'Ecotype')
emm <- emmeans(m1, "Ecotype")
pairs(emm)


### Dcarthusianorum
eschikon.mat<-read.table('01_Dcar_phenoTraits.txt',sep='\t',header=T,na.string='NA')
eschikon.mat$Population <- factor(eschikon.mat$Population, levels = c("Unterstalden","Niedergampel","Grengiols","Gibidumsee","Faldumalp","Simplonpass"))
# log transform
eschikon.mat[,3]<-log(eschikon.mat[,3])
# remove 2SD oultiers
tmp1<-cbind(eschikon.mat[,c(1,3)])
y<-c()
pops<-c("Unterstalden","Niedergampel","Grengiols","Gibidumsee","Faldumalp","Simplonpass") # must be in the same order as in table
for (p in pops){
  tmp2<-subset(tmp1,Population==p)
  x<-tmp2[,2]
  x[x > (mean(x,na.rm=T)+(2*sd(x,na.rm=T)))] <- NA
  x[x < (mean(x,na.rm=T)-(2*sd(x,na.rm=T)))] <- NA
  y<-c(y,x)
}
mat <- data.frame(
  Population = eschikon.mat$Population,
  Flowering_date.d. = y
)
mat$Population <- factor(mat$Population, levels = c("Unterstalden","Niedergampel","Grengiols","Gibidumsee","Faldumalp","Simplonpass"))
# make plots
boxplot(Flowering_date.d. ~ Population, data = mat, col = c('red', 'red', 'red', 'blue', 'blue', 'blue'),las=2)
#Statistics for difference in flowering time between high- and low-elevation
mat<-data.frame(cbind(mat,eschikon.mat$Ecotype,eschikon.mat$Motherline))
colnames(mat)[3:4]<-c('Ecotype','Motherline')
mat$Flowering_date.d.<-as.numeric(as.character(mat$Flowering_date.d.))
m1<-lmer(Flowering_date.d. ~ Ecotype + (1|Population/Motherline), data=mat,  na.action=na.exclude)
m2 <- lmer(Flowering_date.d. ~ 1 + (1|Population/Motherline), data=mat,  na.action=na.exclude)
anova(m1,m2)
lsmeans(m1, 'Ecotype')
emm <- emmeans(m1, "Ecotype")
pairs(emm)




