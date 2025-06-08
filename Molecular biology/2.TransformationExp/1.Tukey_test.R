# This script compares estimates of leaf number and flowering time in the Arabidopsis transformants

library(ggplot2)
library(agricolae)
library(Rmisc)

traits<-c('TLN_SD','days_till_bolting_SD','TLN_LD','days_till_bolting_LD')

mat_mean<-c()
mat_sd<-c()

for (t in traits){
  print(t)
  data<-read.table(paste(t,'.txt',sep=''),header=T,sep='\t',na.string='NA',fill=T)
  df<-data.frame(t(data))
  Geno<-c(rep('High',9),rep('Low',9),rep('Col',2))
  df<-cbind(df,Geno)
  
  High<-subset(df,Geno=='High')[,-9]
  Low<-subset(df,Geno=='Low')[,-9]
  Col<-subset(df,Geno=='Col')[,-9]
  
  vHigh<-as.vector(unlist(High))
  vLow<-as.vector(unlist(Low))
  vCol<-as.vector(unlist(Col))
  x<-c(vHigh,vLow,vCol)
  y<-c(rep('High',length(vHigh)),  rep('Low',length(vLow)),  rep('Col',length(vCol)))
  df2<-data.frame(x,y)
  #boxplot(x~y, data=df2)
  
  m1<-aov(x~y,data=df2)
  m2<-TukeyHSD(m1)
  write.table(m2$y, file=paste(t,'.TukeyHSD.txt',sep=''),quote=F,sep='\t')
  
  results <- HSD.test(m1, "y", group=TRUE)
  results$groups
  
  oder.group <- results$groups[order(rownames(results$groups)),]
  
  stats <- summarySE(df2, measurevar="x", groupvars=c("y"),na.rm=TRUE)
  #stats
  write.table(stats, file=paste(t,'.stats.txt',sep=''),quote=F,row.names=F,sep='\t')
  
  
  p<-ggplot(stats, aes(x=y, y=x)) + 
    geom_bar(position=position_dodge(.5), stat="identity",
             colour="black", width = 0.5) +      # Change bar size with width = 0.5
    geom_errorbar(aes(ymin=x-se, ymax=x+se),
                  size=.5,    # Thinner lines
                  width=.2,
                  position=position_dodge(.5)) + 
    theme_classic() + geom_text(aes(x=y, y=x+se+2,label=oder.group$groups), position=position_dodge(width=0.9), size=4)
  ggsave(paste(t,'.tukey.pdf',sep=''),p) 
}







