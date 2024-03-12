
args <- commandArgs(trailingOnly = TRUE)


data<-read.table(args[1],sep="\t",header=T)
d<-data.frame(data)

pairs<-strsplit(args[2],",")[[1]]
print(pairs)
print(class(pairs))


for (pair in pairs){
	print(paste('state.',pair,sep=''))
	p<-paste('state.',pair,sep='')
	d<- d[ d[[p]] == 1 , ]    
}
#print(head(d))


v<-c()
for (pair in pairs){
	v<-paste(v,pair,sep='_')
}
#print(v)
write.table(d,file=paste(args[3],v,'.txt',sep=''),quote=F,sep='\t',row.name=F)
