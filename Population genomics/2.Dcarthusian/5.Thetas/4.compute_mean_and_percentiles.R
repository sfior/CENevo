pops<-c('pop1','pop2','pop3','pop4','pop5','pop6')
for (p in pops){
	data<-read.table(paste(p,'_w500_s500_filter_250bp.thetasWindow.pestPG',sep=''),header=T,sep='\t')
	qval<-c(.05,.95)
	m_tP<-mean(data[,5])
	q_tP<-quantile(data[,5],qval)
	out_tP <- matrix(data = c(m_tP,q_tP), ncol = 3)
	colnames(out_tP)<-c('mean',qval)
	write.table(out_tP,paste(p,'_mean_and_percentiles.tP',sep=''),sep='\t',row.names=F,col.names=T,quote=F)

	
	m_TajD<-mean(data[,9])
	q_TajD<-quantile(data[,9],qval)
	out_TajD <- matrix(data = c(m_TajD,q_TajD), ncol = 3)
	colnames(out_TajD)<-c('mean',qval)
	write.table(out_TajD,paste(p,'_mean_and_percentiles.TajD',sep=''),sep='\t',row.names=F,col.names=T,quote=F)
}


