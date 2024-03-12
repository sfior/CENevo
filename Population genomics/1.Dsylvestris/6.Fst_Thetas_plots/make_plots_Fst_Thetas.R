# This scripts makes Fst and Theta plots for multiple populations


args <- commandArgs(trailingOnly = TRUE)
if(length(args)<9){
    print("Invalid arguments supplied. Check executuion command.")
    quit()
}



options(warn=1)    
options(scipen=10)	

####    autoloess function    ##################################
####    see https://gist.github.com/kylebgorman/6444612   ######
aicc.loess <- function(fit) {
    # compute AIC_C for a LOESS fit, from:
    # 
    # Hurvich, C.M., Simonoff, J.S., and Tsai, C. L. 1998. Smoothing 
    # parameter selection in nonparametric regression using an improved 
    # Akaike Information Criterion. Journal of the Royal Statistical 
    # Society B 60: 271–293.
    # 
    # @param fit        loess fit
    # @return           'aicc' value
    stopifnot(inherits(fit, 'loess'))
    # parameters
    n <- fit$n
    trace <- fit$trace.hat
    sigma2 <- sum(resid(fit) ^ 2) / (n - 1)
    return(log(sigma2) + 1 + 2 * (2 * (trace + 1)) / (n - trace - 2))
}

autoloess <- function(fit, span=c(as.numeric(args[4]), as.numeric(args[5]))) {
    # compute loess fit which has span minimizes AIC_C
    # 
    # @param fit        loess fit; span parameter value doesn't matter
    # @param span       a two-value vector representing the minimum and 
    #                   maximum span values
    # @return           loess fit with span minimizing the AIC_C function
    stopifnot(inherits(fit, 'loess'), length(span) == 2)
    # loss function in form to be used by optimize
    f <- function(span) aicc.loess(update(fit, span=span))
    # find best loess according to loss function
    return(update(fit, span=optimize(f, span)$minimum))
}
#################################################################

## Fst

if(args[12] == "pdfyes"){
pdf(args[36+11],onefile=T,width = 100, height = 15)
panels<-c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15)
layout(matrix(c(panels,c(15+panels),c(30+panels)),3,15,byrow = TRUE), c(1,1), c(0.6,0.6,0.6,0.6,0.6),TRUE)
par(mar=c(2.1,4.1,0,2.1))
par(omi=c(2,2,2,5))
}

scaff_file<-read.table(args[1],header=F,sep="\t")

for (n in 1:nrow(scaff_file)){
f<-scaff_file[n,][,1]
peakx1<-scaff_file[n,][,2]   #  these are used to make the shaded rectangle of that indicates the Fst peak
peakx2<-scaff_file[n,][,3]
print(f)

raw_data<-read.table(paste('Fst_files/',args[2],"_",f,args[3],sep=""),sep="\t",na.strings=c("na","-0.00000000","0.00000000"))

myvars<-strsplit(args[7],",")[[1]]
Fst<-raw_data[myvars]

pairs<-strsplit(args[8],",")[[1]]
colnames(Fst)<-c("scaffold","base",pairs)

# makes sure that the table is sorted according to ascending base number
Fst<-data.frame(Fst)
Fst<-Fst[with(Fst, order(base)),]


pairs_to_compute<-strsplit(args[9],",")[[1]]
mypairs<-as.list(pairs_to_compute)


for (k in mypairs){
	
	# ---> Creates a 2 column dataframe(i.e. base + Fst values)
	raw_data2<-cbind(Fst[,"base"],Fst[,k])
	data<-na.omit(raw_data2)
	colnames(data)<-c("b","p1")
	data<-data.frame(data)
	#print(head(data))
	b<-data$b
	p1<-data$p1
	#plot(data,cex=0.1,col="grey50")
	
	bnew<-seq(b[1],b[length(b)],as.numeric(args[13]))   # creates a new set of bases to use in the predict function below. Chnange last number to mofify the step of the range
	
	# ---> runs the autoloess to find best fitting 'span' parameter
	set.seed(62485)
	basic <- loess(p1 ~ b)
	good <- autoloess(basic)
	opt_span<-good$pars$span
	
	good_2 <- predict(loess(p1 ~ b,span=opt_span),bnew)  	

	# ---> runs the Permutations
	w<-cbind()
	for (i in 1:as.numeric(args[6])){
		data2<-data.frame(b,sample(p1))
		colnames(data2)<-c("b","p1_2")
		b<-data2$b
		p1_2<-data2$p1_2
		lw1 <- predict(loess(p1_2 ~ b, span=opt_span),bnew)   # this predicts every permutation every N bp
		#lines(data2$b,lw1$fitted,col=2,lwd=2)
		
		w<-cbind(w,lw1)
	}

	wt<-t(w)
	ci<-c()
	for (i in 1:ncol(wt)){
	q<-quantile(wt[,i],c(as.numeric(args[10]),as.numeric(args[11])))
	ci<-rbind(ci,q)
	}
	
	s<-cbind(bnew,ci,good_2)
	
	pops<-unlist(strsplit(k, split='Pair', fixed=TRUE))[2]
	pA<-unlist(strsplit(pops, split='_', fixed=TRUE))[1]
	pB<-unlist(strsplit(pops, split='_', fixed=TRUE))[2]
	meanFst<-read.table(paste('Fst_files/','pop',pA,'.pop',pB,'_mean_and_percentiles.txt',sep=''),sep='\t',header=T)[,1]
	percentileFst<-read.table(paste('Fst_files/','pop',pA,'.pop',pB,'_mean_and_percentiles.txt',sep=''),sep='\t',header=T)[,3]
	#print(percentileFst)
	
	# ---> creates plots on PFD if argument == pdfyes 
	if(args[12] == "pdfyes"){
		

		plot(1,type="n",xlim=c(b[1],b[length(b)]),ylim=c(0,0.8),xlab="",ylab="",axes=F)
		lim <- par("usr")
		rect(peakx1, lim[3]-1, peakx2, lim[4]+1, border = "khaki1", col = "khaki1", density= 30)
		axis(1,labels=F)
		axis(2,labels=T,cex.axis=3)
		box()
		points(b,p1,pch=19,cex=0.5,col="darkgrey")
		points(s[,1],s[,4],type="l",col = 1,lwd=2)
		
		abline(h=meanFst,col='darkgrey',lty=2,lwd=2)
		legend("topleft",paste("span=",round(opt_span,digits=3),sep=""),lwd=2,cex=2,bty='n')

		}
		
	#-----> Build matrix for subsequent work. Contains positions, 5%, 95%, lw1j-lw2j
	write.table(s,file=paste(f,"_",k,"_loess",".txt",sep=""),sep="\t",row.names=F,col.names=T,quote=F)

}
title(main=f,cex.main=4,line=5,outer=T)
mtext(paste('Fst','ThetaD',"Tajima's D",sep='\n\n\n\n\n\n'),line=5,outer=T,cex=3,side=4,las=2)




#############################################################################################################################################

# Pi


pairs_to_compute<-strsplit(args[13+4],",")[[1]]
mypairs<-as.list(pairs_to_compute)
altitudes<-strsplit(args[36+10],",")[[1]]


scaff_file<-read.table(args[13+3],header=F,sep="\t")
scaff_name<-scaff_file[,1]

for (k in mypairs){
				# this is needed to choose colors only:
				print(k)
				index<-match(k,pairs_to_compute)
				alt<-altitudes[index]
				if (alt=='HL'){
					col1=4
					col2=2
				}
				if (alt=='LL'){
					col1=2
					col2="orange"
				}
				if (alt=='HH'){
					col1=4
					col2=5
				}
			if (args[24+11+1]=='Dsyl'){
			if (k=='Pair3_5'){
				k<-'Pair5_3'}
			if (k=='Pair4_5'){
				k<-'Pair5_4'}
				}
			if (args[24+11+1]=='Dcar'){
			if (k=='Pair1_2'){
				k<-'Pair2_1'}
			if (k=='Pair1_4'){
				k<-'Pair4_1'}
			if (k=='Pair1_5'){
				k<-'Pair5_1'}
			if (k=='Pair3_4'){
				k<-'Pair4_3'}
			if (k=='Pair3_5'){
				k<-'Pair5_3'}
				}
			k<-gsub("Pair","",k)
			m = unlist(strsplit(k, split='_', fixed=TRUE))[1]
			n = unlist(strsplit(k, split='_', fixed=TRUE))[2]

			print(f)
			print(paste("pop",m,sep=""))
			print(paste("pop",n,sep=""))
			data1<-read.table(paste('tP_files/',args[13+11],"pop",m,args[13+1],"_",f,args[13+2],sep=""),sep="\t", header=F,na.string="0.00000000")
			data2<-read.table(paste('tP_files/',args[13+11],"pop",n,args[13+1],"_",f,args[13+2],sep=""),sep="\t", header=F,na.string="0.00000000")



o<-merge(data1,data2, by.x = "V2", by.y="V2", all.y = TRUE, all.x=TRUE)

#  ------------------------> permute datasets, make 2 lowess, substract

b<-o[,1]
p1<-o[,3]
p2<-o[,6]

bj<-b
#----> subtraction of the lines with loess
data<-data.frame(b,p1)
set.seed(62485)
basic <- loess(p1 ~ b)
good <- autoloess(basic)
opt_span1<-good$pars$span

good_2 <- predict(loess(p1 ~ b,span=opt_span1),data.frame(bases = b)) 
	lw1j<-good_2



data<-data.frame(b,p2)
set.seed(62485)
basic <- loess(p2 ~ b)
good <- autoloess(basic)
opt_span2<-good$pars$span

good_2 <- predict(loess(p2 ~ b,span=opt_span2),data.frame(bases = b)) 
	lw2j<-good_2


pA<-unlist(strsplit(k, split='_', fixed=TRUE))[1]
pB<-unlist(strsplit(k, split='_', fixed=TRUE))[2]
meanPi_pA<-read.table(paste('tP_files/','pop',pA,'_mean_and_percentiles.tP',sep=''),sep='\t',header=T)[,1]
meanPi_pB<-read.table(paste('tP_files/','pop',pB,'_mean_and_percentiles.tP',sep=''),sep='\t',header=T)[,1]


plot(1,type="n",xlim=c(bj[1],bj[length(bj)]),ylim=c(-6,15),cex.axis=3,xlab="",ylab="",axes=F)
lim <- par("usr")
rect(peakx1, lim[3]-1, peakx2, lim[4]+1, border = "khaki1", col = "khaki1", density= 30)


abline(h=meanPi_pA,col=col1,lty=2,lwd=0.5)

axis(1,labels=F)
ylab=seq(-9, 22)
axis(2,at=ylab,labels=ylab,cex.axis=3)
axis(4,at=ylab,labels=ylab+2,cex.axis=3,mgp = c(2, 2, 0))
box() 
abline(h=meanPi_pB,col=col2,lty=2,lwd=0.5)

lines(bj,lw1j,col=col1,lwd=2)
lines(bj,lw2j,col=col2,lwd=2)
legend('topleft',c(paste("span=",round(opt_span1,digits=3)),paste("span=",round(opt_span2,digits=3))),lwd=2,cex=2,col=c(1,2),bty='n')


# --------- DIFF --------
d<-log10(lw1j/lw2j)
d<-(d*10)-2

r1<-quantile(d,.95,na.rm=T)
r2<-quantile(d,.05,na.rm=T)
z<-cbind(bj,d)

#-------> Permutation of original datasets and subtraction of loess
w<-cbind()
for (i in 1:args[13+6]){
perm1<-data.frame(bj,sample(p1))
colnames(perm1)<-c("bj","p1s")
perm2<-data.frame(bj,sample(p2))
colnames(perm2)<-c("bj","p2s")

lw1 <- predict(loess(perm1$p1s ~ b, span=opt_span1), data.frame(bases = b))
lw2 <- predict(loess(perm2$p2s ~ b, span=opt_span2), data.frame(bases = b))
lw1j<-lw1
lw2j<-lw2


# --------- DIFF --------
ds<-log10(lw1j/lw2j)
ds<-(ds*10)-2
w<-cbind(w,ds)
}


# ------> confidence intervals
wt<-t(w)
ci<-c()
for (i in 1:ncol(wt)){
q<-quantile(wt[,i],c(as.numeric(args[13+7]),as.numeric(args[13+8])),na.rm=TRUE)
ci<-rbind(ci,q)
}

#-----> Build matrix for subsequent work. Contains positions, 5%, 95%, lw1j-lw2j
s<-cbind(bj,ci,d)
write.table(s,file=paste(f,"_Pop_",m,"_",n,args[13+2],".txt",sep=""),sep="\t",row.names=F,col.names=T,quote=F)

colnames(s)<-c("bj","ci1","ci2","d")
s<-data.frame(s)
lines(d~bj,data=s,col='darkgreen',lwd=2)
lines(ci1~bj,data=s,col='darkgrey',lwd=1)
lines(ci2~bj,data=s,col='darkgrey',lwd=1)

		}  
	

#######################################################################################################################

##  Tajima

pairs_to_compute<-strsplit(args[24+4],",")[[1]]
mypairs<-as.list(pairs_to_compute)
altitudes<-strsplit(args[36+10],",")[[1]]


for (k in mypairs){
				# this is needed to choose colors only:
				print(k)
				index<-match(k,pairs_to_compute)
				alt<-altitudes[index]
				if (alt=='HL'){
					col1=4
					col2=2
				}
				if (alt=='LL'){
					col1=2
					col2="orange"
				}
				if (alt=='HH'){
					col1=4
					col2=5
				}

			if (args[24+11+1]=='Dsyl'){
			if (k=='Pair3_5'){
				k<-'Pair5_3'}
			if (k=='Pair4_5'){
				k<-'Pair5_4'}
				}
			if (args[24+11+1]=='Dcar'){
			if (k=='Pair1_2'){
				k<-'Pair2_1'}
			if (k=='Pair1_4'){
				k<-'Pair4_1'}
			if (k=='Pair1_5'){
				k<-'Pair5_1'}
			if (k=='Pair3_4'){
				k<-'Pair4_3'}
			if (k=='Pair3_5'){
				k<-'Pair5_3'}
				}
			k<-gsub("Pair","",k)
			m = unlist(strsplit(k, split='_', fixed=TRUE))[1]
			n = unlist(strsplit(k, split='_', fixed=TRUE))[2]

			print(f)
			print(paste("pop",m,sep=""))
			print(paste("pop",n,sep=""))
			data1<-read.table(paste('Tajima_files/',args[24+11],"pop",m,args[24+1],"_",f,args[24+2],sep=""),sep="\t", header=F,na.string="0.00000000")
			data2<-read.table(paste('Tajima_files/',args[24+11],"pop",n,args[24+1],"_",f,args[24+2],sep=""),sep="\t", header=F,na.string="0.00000000")

o<-merge(data1,data2, by.x = "V2", by.y="V2", all.y = TRUE, all.x=TRUE)

#  ------------------------>  permute datasets, make 2 lowess each time, substract

b<-o[,1]
p1<-o[,3]
p2<-o[,6]

bj<-b
#----> subtraction of the lines with loess
data<-data.frame(b,p1)
set.seed(62485)
basic <- loess(p1 ~ b)
good <- autoloess(basic)
opt_span1<-good$pars$span

good_2 <- predict(loess(p1 ~ b,span=opt_span1),data.frame(bases = b)) 
	lw1j<-good_2



data<-data.frame(b,p2)
set.seed(62485)
basic <- loess(p2 ~ b)
good <- autoloess(basic)
opt_span2<-good$pars$span

good_2 <- predict(loess(p2 ~ b,span=opt_span2),data.frame(bases = b)) 
	lw2j<-good_2

pA<-unlist(strsplit(k, split='_', fixed=TRUE))[1]
pB<-unlist(strsplit(k, split='_', fixed=TRUE))[2]
meanTajD_pA<-read.table(paste('Tajima_files/','pop',pA,'_mean_and_percentiles.TajD',sep=''),sep='\t',header=T)[,1]
meanTajD_pB<-read.table(paste('Tajima_files/','pop',pB,'_mean_and_percentiles.TajD',sep=''),sep='\t',header=T)[,1]


plot(1,type="n",xlim=c(bj[1],bj[length(bj)]),ylim=c(-7,2),cex.axis=3,xlab="",ylab="",axes=F)
lim <- par("usr")
rect(peakx1, lim[3]-1, peakx2, lim[4]+1, border = "khaki1", col = "khaki1", density= 30)

abline(h=meanTajD_pA,col=col1,lty=2,lwd=0.5)


axis(1,labels=F)
ylab=seq(-9,4)
axis(1,labels=T,cex.axis=3,mgp = c(2, 2, 0))
axis(2,at=ylab,labels=ylab,cex.axis=3)
axis(4,at=ylab,labels=ylab+5,cex.axis=3,mgp = c(2, 2, 0))
box() 
mtext(paste("Pop",m,"vsPop",n,sep=''), side = 1, line = 10, cex=2)


abline(h=meanTajD_pB,col=col2,lty=2,lwd=0.5)

lines(bj,lw1j,col=col1,lwd=2)
lines(bj,lw2j,col=col2,lwd=2)
legend('topleft',c(paste("span=",round(opt_span1,digits=3)),paste("span=",round(opt_span2,digits=3))),lwd=2,cex=2,col=c(1,2),bty='n')


# --------- DIFF --------

d<-lw1j-lw2j
d<-(d*2)-5


r1<-quantile(d,.95,na.rm=T)
r2<-quantile(d,.05,na.rm=T)
z<-cbind(bj,d)


#-------> Permutation of original datasets and subtraction of loess
w<-cbind()
for (i in 1:args[24+6]){
perm1<-data.frame(bj,sample(p1))
colnames(perm1)<-c("bj","p1s")
perm2<-data.frame(bj,sample(p2))
colnames(perm2)<-c("bj","p2s")


lw1 <- predict(loess(perm1$p1s ~ b, span=opt_span1), data.frame(bases = b))
lw2 <- predict(loess(perm2$p2s ~ b, span=opt_span2), data.frame(bases = b))
lw1j<-lw1
lw2j<-lw2



# --------- DIFF --------

ds<-lw1j-lw2j
ds<-(ds*2)-5


w<-cbind(w,ds)
}


# ------> confidence intervals
wt<-t(w)
ci<-c()
for (i in 1:ncol(wt)){
q<-quantile(wt[,i],c(as.numeric(args[24+7]),as.numeric(args[24+8])),na.rm=TRUE)
ci<-rbind(ci,q)
}

#-----> Build matrix for subsequent work. Contains positions, 5%, 95%, lw1j-lw2j
s<-cbind(bj,ci,d)
write.table(s,file=paste(f,"_Pop_",m,"_",n,args[24+2],".txt",sep=""),sep="\t",row.names=F,col.names=T,quote=F)

colnames(s)<-c("bj","ci1","ci2","d")
s<-data.frame(s)
lines(d~bj,data=s,col = 'darkgreen',lwd=2)
lines(ci1~bj,data=s,col='darkgrey',lwd=1)
lines(ci2~bj,data=s,col='darkgrey',lwd=1)

		}  
	

	
	}




if(args[12] == "pdfyes"){
dev.off()
}	