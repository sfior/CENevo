# This script combines the fst values of each pair of pops in a single table. 
# It is written to combine multiple pairs produced by angsd

# Syntax from bash shell:
# Rscript combine_fst_from_multiple_angsd_pairs.R \
# FILE_WITH_LIST_OF_PAIRS \
# SLIDING WINDOW SIZE \
# STEP SIZE
# PERCENTAGE OF COVERAGE FOR WINDOW TO BE RETAINED

# Example from bash shell:
# Rscript combine_fst_from_multiple_angsd_pairs.R fst_comps.txt 500 500 50

args <- commandArgs(trailingOnly = TRUE)

sw<-args[2]
step<-args[3]
ext=paste('_w',sw,'_s',step,'.fst',sep='')

COMPS<-read.table(args[1])
COMPS<-COMPS[,1]
print(COMPS)


for (comp in COMPS){
	data<-read.table(paste(comp,ext,sep=''),header=T,sep='\t',row.names=NULL)
	colnames(data)<-c('region','chr','midPos','Nsites','fst')
#	print(data)
	if (exists('m1')){
		m1<-merge(m1,data,by=c('region','chr','midPos','Nsites'))
		} else {
		m1<-data
		}
	}	
colnames(m1)<-c('region','chr','midPos','Nsites',paste(COMPS,sep=','))
#print(m1)
write.table(m1,file='combined_pairs_unfiltered.fst',row.names = F,col.names = T,quote=F,sep='\t')

sw<-as.numeric(args[2])
m1<-data.frame(m1)
if (sw > 1){
filter<-(sw/100)*as.numeric(args[4])
print(paste('********** Applying filter > ',filter,' ************',sep=''))
m2<-subset(m1,Nsites > filter)
#print(m2)
}
write.table(m2,file=paste('combined_pairs_filter_',filter,'bp.fst',sep=''),row.names = F,col.names = T,quote=F,sep='\t')

m3<-m2[-1]
m3<-m3[with(m3, order(chr,midPos)),]
#print(m3)

write.table(m3,file=paste('combined_pairs_filter_',filter,'bp_sorted.fst',sep=''),row.names = F,col.names = T,quote=F,sep='\t')


warnings()



