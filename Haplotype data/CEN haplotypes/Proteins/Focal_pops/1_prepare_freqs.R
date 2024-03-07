# This scripts calculates frequencies of protein variants for each population 

coords<-read.table('coordinates.txt',header=F,sep='\t')
data<-read.table('Proteins_alignm_converted_ATCG.phy',header=F,sep=' ',fill=T)[-1,][,-3]

pops<-coords[,1]
uniseqs<-sort(unique(data[,2]))

new<-cbind(coords,matrix(nrow=nrow(coords),ncol=length(uniseqs)+1))

for (n in 1:nrow(new)){
	p<-as.character(new[n,][1][[1]])
	q<-length(grep(p,data[,1]))
	new[n,][4]<-q
}
print(new)

for (m in 1:length(uniseqs)){
	seq=uniseqs[m]
	for (n in 1:nrow(new)){
		p<-as.character(coords[n,][1][[1]])
		tmp<-grep(p,data[,1])
		#print(data[tmp,])
		tmp2<-data[tmp,]
		#print(tmp2)
		
		q<-length(grep(seq,tmp2[,2]))
		print(q)
		new[n,][4+m]<-q/new[n,][4]
	} 
}

colnames(new)<-c('pop','Latitude','Longitude','Tot_haplos','freq_h1','freq_h2','freq_h3','freq_h4','freq_h5','freq_h6')
#write.table(new,file='Proteins_freqs.txt',row.names=F,quote=F,sep='\t')




