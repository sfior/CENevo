## This script computes 1D SFS from the MLE of the Allele Freqs as output by poolstat (in the _sites  file)

## Syntax
#Rscript make_1DSFS_poolstat.R   poolstatoutput.1-20_sites_lowerQ20eq0.txt

args <- commandArgs(trailingOnly = TRUE)

nnorm <- function(x) x/sum(x)


data<-read.table(args[1],header=F,sep="\t")
data2<-cbind(data[,c(6:11)])
ac<-data2*40
ac<-data.frame(ac)
colnames(ac)<-c(1:6)

#### make 1D SFS
sfs<-c()
for (p in 1:6){
	pop.sfs<-c()
	mat<-ac[p]
	#print(head(mat))
	colnames(mat)<-c('X')
	for (n in 0:40){ 
		print(n)
		v<-subset(mat,X==n)
		print(nrow(v))
		pop.sfs<-c(pop.sfs,nrow(v))
		}
		print(pop.sfs)
		
		#write sfs to files 
		cat(pop.sfs,file=paste('pop',p,'.sfs',sep=''))
		res2<-rbind(pop.sfs)
		res2 <- t(apply(res2,1,nnorm))
		cat(res2,file=paste('pop',p,'.norm.sfs',sep=''))
		
		# add pop.sfs to 6 pop sfs
		sfs<-rbind(sfs,pop.sfs)		
	}
print(sfs)
pdf(file='1D.SFS.6pops.pdf')
res<-sfs[,-1]
colnames(res) <- 1:40
# density instead of expected counts
res <- t(apply(res,1,nnorm))
#plot the none ancestral sites
barplot(res,beside=T,legend=c("POP1","POP2","POP3","POP4","POP5","POP6"),names=1:40,main="realSFS non ancestral sites")
dev.off()


