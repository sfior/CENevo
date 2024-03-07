# This scipt plots frequencies of protein variants as pie charts on a map


library(maps)
library(mapdata)
library(maptools)  #for shapefiles
library(scales)  #for transparency
library(mapplots)

data<-read.table("Protein_freqs.txt",header=T,sep='\t')

# Move some charts to make them more visible
x<-subset(data,Population=='Lac_de_Mauvoisin')
x[2]<-47.8  #Lac
x[3]<-6.5
data[row.names(x),]<-x

x<-subset(data,Population=='Challer')
x[2]<-47.8  #Challer
x[3]<-7
data[row.names(x),]<-x

x<-subset(data,Population=='Tsanfleuron')
x[2]<-47.8  #Tsanfleuron
x[3]<-7.5
data[row.names(x),]<-x

x<-subset(data,Population=='Saxon')
x[2]<-47.8  #Saxon
x[3]<-8
data[row.names(x),]<-x

x<-subset(data,Population=='Saviese')
x[2]<-47.8  #Saviese
x[3]<-8.5
data[row.names(x),]<-x

x<-subset(data,Population=='Varen')
x[2]<-47.8  #Varens
x[3]<-9
data[row.names(x),]<-x

x<-subset(data,Population=='Rifugio_Graziani')
x[2]<-45.5 
x[3]<-10
data[row.names(x),]<-x

x<-subset(data,Population=='Forte_San_Marco_Low')
x[2]<-45.5 
x[3]<-10.5
data[row.names(x),]<-x

x<-subset(data,Population=='Forte_San_Marco_High')
x[2]<-46
x[3]<-11
data[row.names(x),]<-x

x<-subset(data,Population=='Polsa')
x[2]<-45.7 
x[3]<-11.3
data[row.names(x),]<-x

x<-subset(data,Population=='Soca_High')
x[2]<-46.2 
x[3]<-13.3
data[row.names(x),]<-x

x<-subset(data,Population=='Soca_Low')
x[2]<-46.2 
x[3]<-13.7
data[row.names(x),]<-x



dev.new(width=15, height=15)
#pdf(file='Rplot.pdf',width=15, height=15)
my_col<-c(alpha('blue', 0.6), alpha('red',0.6))

map("world", xlim=c(5,20),ylim=c(38,48), col="gray90", fill=TRUE) 
for (n in 1:nrow(data)){
	#print(n)
	line<-data[n,]
	add.pie(z=c(line[,8],line[,9]), x=line[,3], y=as.numeric(line[,2]), radius=sqrt(line[,5]/500), col=my_col, labels='')
}
#dev.off()

# ==========================
# Test correlations in Mte Baldo

Baldo<-subset(data, Population=='Lumini' |  Population=='Saviese' | Population=='Lac_de_Mauvoisin' | Population=='Forte_San_Marco_High' | Population=='Rifugio_Graziani' | Population=='Saxon' | Population=='Challer' | Population=='Tsanfleuron' | Population=='Varen' | Population=='Polsa' | Population=='Forte_San_Marco_Low' )
p<-as.character(Baldo$Popoulation)
x<-as.numeric(Baldo$Altitude)
y<-as.numeric(Baldo$FreqHigh)
plot(x,y)
with(Baldo, text(FreqHigh~Altitude, labels = Baldo$Population, pos = 1))
abline(fit<-lm(y~x))
corr<-cor.test(x,y, method='spearman',alternative = "g",exact = FALSE)
legend("topleft", bty="n", legend=paste('y =', format(coef(fit)[[2]],digits=4), 'x', '+', format(coef(fit)[[1]],digits=4),"\n","R2 = ",format(summary(fit)$adj.r.squared, digits=4),"\n","Spearman rho = ",format(corr[4]),"\n","pvalue = ",format(corr[3]),sep=""))




