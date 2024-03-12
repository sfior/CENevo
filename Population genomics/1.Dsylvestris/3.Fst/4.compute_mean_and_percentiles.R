data<-read.table('combined_pairs_filter_250bp_sorted.fst',header=T,sep='\t')
pairs<-read.table('fst_comps.txt',header=F)
pairs<-pairs[,1]

for (p in pairs){
	d<-data[[p]]
	m<-(mean(d))
	qval<-c(.05,.9,.95)
	q<-quantile(d,qval)
	out <- matrix(data = c(m,q), ncol = 4)
	colnames(out)<-c('mean',qval)
	write.table(out,paste(p,'_mean_and_percentiles.txt',sep=''),sep='\t',row.names=F,col.names=T,quote=F)
	}

