# This script performs partial Mantel test on a genome scan to test environmental associations, using a third matrix to correct for structure
# It takes allele freqs for every SNP position in a scaffold and performs env ass with a a number of variables provided in a table

# It takes the following arguments:
# 1: input file for genetic data (poosltat output subset for the target scaffolds)
# 2: name of the scaffold analysed, only for output file
# 3: input file for environmental data (generally with multiple variables)
# 4: the order in which the populations of the envritomental data must be rearranged in orde to be IDENTICAL to the order 
#		of the pops in the genetic data
#		Eg: 6,4,2,5,1,3
# 5: list of column names for the genetic data. Scaffold must be named 'chr' and SNP position must be names 'pos' (used for subsetting)
#		Eg:  "chr","pos","Tsanfleuron_POP1.txt","Challer_POP2.txt","Saviese_POP3.txt","Saxon_POP4.txt","LacdeMauv_POP5.txt","Varen_POP6.txt"
# 6: wheater to subset the gentic table. If so, it must be 'yes' 
# 7: minimum base postion for subsetting (e.g. 200000)  
# 8: maxumum base position for subsetting (e.g. 220000)   
# 9: column number of environmental data to analyse. Can be a value or passed from IDX
# 10: number of replicates in mantel test 
# 11: File name of third matrix


# Examples: 
# With subsetting only a portion of the scaffold:
#Rscript Mantel_test_alleles_variables_v3.R \
#scaffold1_size1318325_1_614794 \
#scaffold1 \
#Dsylvestris_all.txt \
#6,4,2,5,1,3 \
#"chr","pos","Tsanfleuron_POP1.txt","Challer_POP2.txt","Saviese_POP3.txt","Saxon_POP4.txt","LacdeMauv_POP5.txt","Varen_POP6.txt" \
#yes 700 800    ->  In case you don't want to subset, write a large span (e.g 1 1000000)
# 6 999
# global.fst.matrix.txt




library(vegan)

args <- commandArgs(trailingOnly = TRUE)


data <- read.table(args[1], header=F,sep="\t",na.string="NA")
scaff<-args[2]
var <- read.table(args[3], header=T,sep="\t")
sort<-strsplit(args[4],",")[[1]]
sort<-rev(sort)
var_sorted<-c()
for (so in sort){
	var_sorted<-rbind(var[so,],var_sorted)	
}
write.table(var_sorted,file="Environmental_variables_ordered.txt",sep="\t",col.names=T,row.names=F,quote=F)

Fst_matrix<-read.table(args[11], header=T,sep="\t")
Fst.dist<-as.dist(Fst_matrix[,-1])


cols<-strsplit(args[5],",")[[1]]
colnames(data)<-c(cols)
d<-data.frame(data)
d<-na.omit(d)

if(args[6]=='yes'){  
	min<-as.numeric(args[7])
	max<-as.numeric(args[8])
    outliers<-subset(d,pos > min & pos < max )
	d<-outliers
}

print('##########   start associations  ############')

n=as.numeric(args[9])
	print(colnames(var_sorted)[n])
	out<-c()
	for (s in 1:nrow(d)){
		s<-d[s,]  ## change SNP
		print(s[1:2])
		f<-s[3:8]
		if (var(as.numeric(f[1,])) != 0) {   # this is needed for SNPs that have equal AF across all pops, as this breaks the mantel test. It tests if all values are NOT the same.
		
		ft<-t(f)
		v<-var_sorted[n]   ## change var
		m<-cbind(ft,v)
		#print(m)
		af.vegdist<-vegdist(m[1],method='euclidean')
		var.vegdist<-vegdist(m[2],method='euclidean')
		part_mtest<-mantel.partial(af.vegdist,var.vegdist,Fst.dist,method='pearson',permutations=as.numeric(args[10]))
		l<-cbind(s,part_mtest$signif)
		out<-rbind(l,out)
		}
	}	
	write.table(out,file=paste(scaff,'_',colnames(var_sorted)[n],'.txt',sep=""),sep="\t",col.names=T,row.names=F,quote=F)




