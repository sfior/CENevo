# This scipt makes a PCA of the focal populations in environmental space 

library(ggbiplot)

vars<-c("mind","prec","srad","tave","clou","ept","mbal","tmax","tmin","cdnx","ddeg000","ddeg300","ddeg556","exposi","alti","slope","pday","sfroyy","swb","topos","twi25ss")

data2<-read.table('Dcarthusian_all.txt',header=T,sep='\t')
data_env2<-data2[6:26]
colnames(data_env2)<-vars
env.pca2 <- prcomp(na.omit(data_env2), center = T, scale. = TRUE) 
g <- ggbiplot(env.pca2, obs.scale = 1, var.scale = 1, ellipse = F, circle = T, labels =  data2[,1], labels.size = 4) 
g <- g + scale_color_discrete(name = '') + theme(legend.direction = 'horizontal', legend.position = 'top') + theme_bw() 
print(g)
#ggsave("Dcar_ggplot.pdf")





