
# This script performs an HMM analyses to detect regions of high differentiations within scaffolds. 
# The analyses detects swtitches between:
# state1: low_differentiation state -> normal distribution with mean and sd correspoding to those calculate from the actual genome-wide data 
# state2: high_differentiation state -> normal distribution with mean at a certain quantile (defined by the user), and st dev corresponding to the genome-wide value (i.e. equals above)

# Key references:
# Riesch... Nosil 2017: 'Transitions between phases of genomic differentiation during stick-insect speciation'. Nat. Ecol. Evol.
# Hofer... Excoffier 2012: ' Evolutionary forces shaping genomic islands of population differentiation in humans' BMC Genomics


# REQUIRED ARGUMENTS:
#arg1: combined table of fst values
#arg2: column name (of combined table) on which to run analysis (e.g. pop pairs for fst)
#arg3: quantile of the high-diff state distribition (e.g. .9)
#arg4: file with scaffold list on which to compute the BaumWelch. Typically contains either all scaffolds or list of scaffolds in a chromosome


# EXAMPLE
# Rscript  3_compute.HMM.normal.twostates.Mstep.normc.v2.R combined_pairs_filter_250bp_sorted.fst "pop1.pop2","pop1.pop3" .9   chr_13


library("HiddenMarkov")
Mstep.normc <- function(x, cond, pm, pn){
nms <- sort(names(pm))
         n <- length(x)
         m <- ncol(cond$u)
         if (all(nms==c("mean", "sd"))){
             mean <- pm$mean
             sd <- pm$sd
             return(list(mean=mean, sd=sd))
         }
         if (all(nms == "mean")) {
             mean <- pm$mean
             return(list(mean = mean))
         }
         if (all(nms == "sd")) {
             sd <- pm$sd
             return(list(sd = sd))
         }
         stop("Invalid specification of parameters")
     }
rnormc <- rnorm
dnormc <- dnorm
pnormc <- pnorm
qnormc <- qnorm



args <- commandArgs(trailingOnly = TRUE)
if(length(args)<4){
    print("Invalid arguments supplied. Check executuion command.")
    quit()
}

pops<-strsplit(args[2],",")[[1]]
print(pops)

SCAFFSbaum<-scan(args[4],what=character())

SCAFFSviterbi<-scan(args[6],what=character())

all_data<-read.table(args[1],sep="\t",header=T)
all_data<-data.frame(all_data)
print(head(all_data))

# creates the dataset choosing scaffolds for the BaumWlech, maintaining order of scaffolds on chromosomes
data<-c()
for (s in SCAFFSbaum){
	print(s)
	d<-subset(all_data,chr==s)
	d<-d[order(d$midPos),]     # makes sure that positions within scaffolds are in increasing order
	data<-rbind(data,d)
	#print(tail(data))
}
print(head(data))

#------------------------------------------------------------------------------------
# step1: I run dthmm using the distribution of fst from specified scaffolds
#		 state1 has a normal distribution with mean and sd as inferred from the actual data
# 		 state2 has mean at the 90 percentile of the fst dist.  standard dev = the standard dev of the actual fst values
# 		 On this I estimate the dmm parameters using BaumWelch algorithm
#-------------------------------------------------------------------------------------


# Transforms fst data in logit
for (pop in pops){
print(pop)
fst<-data[[pop]]
posfst<-fst[fst > 0]
q<-quantile(posfst,c(.01)) 
fst[fst < q] <- q
fst<-qlogis(fst)

M<-mean(fst)
SD<-sd(fst)
T<-quantile(fst,as.numeric(args[3]))

# Set the m X m transition matrix (Pi). These values will also be estimated from the data
Pi <- matrix(c(0.9, 0.1, 0.1, 0.9), byrow=TRUE, nrow=2)  
delta <- c(0, 1) 


# Run dthmm 
x<-dthmm(fst, Pi, delta, "normc", list(mean=c(T, M), sd=c(SD, SD)),discrete=FALSE)

# run BaumWelch to estimate hmm parameters
y<-BaumWelch(x,bwcontrol(maxiter = 500, tol = 1e-04, prt = F, posdiff = FALSE,converge = expression(diff < tol)))


##### get some stats	
#print(y$Pi)
print(summary(y))
print(logLik(y))
#hist(residuals(y))
##### check parameter estimates
#print(sum(y$delta))
#print(y$Pi %*% rep(1, ncol(y$Pi)))


#--------------------------------------------------------------------------------------
# step2: I use the Pi estimated above to run dthmm
# Viterbi algorithm to infer switches between states
# Becasue dtmm detects switches among consecutive values, this analysis is performed independently within each scaffold 
#--------------------------------------------------------------------------------------

mat<-c()
	

for (s in SCAFFSviterbi){
	print(s)
	scaff_mat<-c()
	ds<-subset(data,chr==s)
	ds<-ds[order(ds$midPos),]    # makes sure that positions within scaffolds are in increasing order
	
	if (nrow(ds) > 20){    # I set the minimum N of windows to 20, which means that I consider only scaffolds at leat 10'000 bp long (20*500=10'000)
	scaff_fst<-ds[[pop]]
	
	scaff_fst[scaff_fst < q] <- q
	scaff_fst<-qlogis(scaff_fst)
	
	x <- dthmm(scaff_fst,  y$Pi, delta, "normc", list(mean=c(T, M), sd=c(SD, SD)),discrete=FALSE)
	#print(summary(x))

	states <- Viterbi(x)
	
	efst<-1/(1+exp(-1 * scaff_fst)) 
	
	scaff_mat<-cbind(ds[,c(1,2,3)],ds[[pop]])
	scaff_mat<-cbind(scaff_mat,efst)
	scaff_mat<-cbind(scaff_mat,states)
	
	mat<-rbind(mat,scaff_mat)
}
}



colnames(mat)<-c('chr','midPos','Nsites',pop,paste('efst.',pop,sep=''),paste('state.',pop,sep=''))
write.table(mat,file=paste(pop,"_",args[4],'.hmm.txt',sep=''),quote=F,sep="\t",row.names=F)	

# a table that combines:
#row1: the re-transformed mean and 90% of the fst dist
#row2: the re-transformed sd of the fst dist
#row3: the proportion of states1 and 2  -> not informative in case of by-scaffold, as refers only to last scaffold of loop
#row4 + row5: the transition matrix of y
out<-round(rbind(1/(1+ exp(-1*y$pm$mean)),1/(1+ exp(-1*y$pm$sd)),table(states)/length(states),y$Pi),4)
write.table(out,file=paste(pop,"_",args[4],".Summarypar.txt",sep=""),row.names=F,col.names=F)

