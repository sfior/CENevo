
## This script computes the2D SFS starting from the MLE of the Allele Freqs as output by poolstat (in the _sites  file)

## Syntax
#Rscipt make_SFS_plot_poolstat_sites.R column1_in_poolstat_file column1_in_poolstat_file

args <- commandArgs(trailingOnly = TRUE)

nnorm <- function(x) x/sum(x)


data<-read.table(args[1],header=F,sep="\t")
data2<-cbind(data[,c(6:11)])
ac<-data2*40
ac<-data.frame(ac)
colnames(ac)<-c(1:6)

#### make 2D SFS

w<-args[2]
y<-args[3]

		mat<-ac[,c(w,y)]
		colnames(mat)<-c('W','Y')
		#print(head(mat))
		twod.sfs<-c()
		for (a in 0:40){
			for (b in 0:40){
				print(paste(a,b))
				bin1<-subset(mat,W==a)
				
				bin2<-subset(bin1,Y==b)
				#print(head(bin2))
				#print(nrow(bin2))
				twod.sfs<-c(twod.sfs,nrow(bin2))
				#print('pippo')
				#print(twod.sfs)
			}
			
		}
		#print(sfs)
		#print(length(sfs))
		
		#write sfs to files 
		cat(twod.sfs,file=paste('pop',w,'.pop',y,'.sfs',sep=''))
		
		res2<-rbind(twod.sfs)
		res2 <- t(apply(res2,1,nnorm))
		cat(res2,file=paste('pop',w,'.pop',y,'.norm.sfs',sep=''))
		

