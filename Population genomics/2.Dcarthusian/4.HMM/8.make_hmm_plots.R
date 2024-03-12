
# Thi script produces either a PNG file with combined manhattan plots for 15 PW, or single PDFs for populations (as specified from the command line)

  

# 1: 'Dsyl' or 'Dcar'
# 2: Pops for the 9 PW comps between altitude
# 3: Pops for the 6 PW comps within altitude 
# 4: png for 15 PW, or pdf for a single pair.  In the latter case, still write something for #3, although not used, else the script breaks

args <- commandArgs(trailingOnly = TRUE)


args <- commandArgs(trailingOnly = TRUE)
if(length(args)<4){
    print("Invalid arguments supplied. Check executuion command.")
    quit()
}


## get some parameters from the command line
sp<-args[1]
pops1<-strsplit(args[2],",")[[1]]
print(pops1)
pops2<-strsplit(args[3],",")[[1]]
print(pops2)

## create the gap between scaffolds to plot on a chromosome
gapsize<-1000
empty<-c(rep(NA,6))
gap<-c()
for (n in 1:gapsize){
gap<-rbind(gap,empty)
}



if (args[4] == 'png'){
png(file=paste("combined_plots_Fst",".png",sep=""),width=5000,height=6000)
layout(matrix(c(1:30),15,2,byrow = TRUE), widths=c(10,1.5),heights=c(1,1),TRUE) 
par(mar=c(12,15,1,1))
par(oma=c(20,10,10,10))


# ------------------------------------------------------------------------------------------------------
#  Here I plot the 1rst PW across altitude to write candidate numbers on top
# ------------------------------------------------------------------------------------------------------
for (p in pops1[1]){
print(p)

mat<-c()
for (n in seq(1:14)){
	data<-read.table(paste('hmm.chr',n,'/',p,'_chr_',n,'_3rd_state.hmm.txt',sep=''),sep='\t',header=T)
	colnames(gap)<-colnames(data)
	mat<-rbind(mat,data,gap)	
}
n<-15
data<-read.table(paste('hmm.chr',n,'/',p,'_chr_',n,'_3rd_state.hmm.txt',sep=''),sep='\t',header=T)
mat<-rbind(mat,data)

p1<-unlist(strsplit(p,".",fixed=T))[1]
p1<-unlist(strsplit(p1,"pop",fixed=T))[2]
p2<-unlist(strsplit(p,".",fixed=T))[2]
p2<-unlist(strsplit(p2,"pop",fixed=T))[2]
mycolors<-c('orange',"darkgrey",'red')
mypch<-c(1,1,19)
mycex<-c(3,3,4)
if (args[1]=='Dsyl'){
plot(1,type="n",xlim=c(1,length(mat[,2])),ylim=c(0,0.8),xlab="",ylab=paste('Pop',p1,'.Pop',p2,sep=''),cex.lab=6,axes=F,mgp=c(10,4,0))
mtext(c(1:13),at=c(6500,9500,30500,36000,38000,39200,75500,94000,117500,123500,131500,136500,142000),cex=3,side=3,col='red')
}
if (args[1]=='Dcar'){
plot(1,type="n",xlim=c(1,length(mat[,2])),ylim=c(0,1),xlab="",ylab=paste('Pop',p1,'.Pop',p2,sep=''),cex.lab=6,axes=F,mgp=c(10,4,0))
mtext(c(1:6),at=c(67000,73000,75500,77000,90500,153000),cex=3,side=3,col='red')
}


points(1:length(mat[,2]),mat[,5],col=mycolors[mat[,6]],pch=mypch[mat[,6]],cex=mycex[mat[,6]])
axis(1,labels=F,tck=0,lwd=6)
axis(2,labels=T,cex.axis=5,pos=-2000,mgp=c(3,2,0),lwd=6)

x<-paste('efst.',p,sep='')
mat2<-subset(mat,get(x)>0.005)
hist(mat2[,5],breaks=1000,main='',xlab="",ylab='',xlim=c(0,0.6),axes=F,col='black')
axis(1,labels=F,lwd=6)
axis(2,labels=T,cex.axis=5,lwd=6)

}



# ------------------------------------------------------------------------------------------------------
#  Here I plot the other 8 PW across altitude
# ------------------------------------------------------------------------------------------------------
for (p in pops1[-1]){
print(p)

mat<-c()
for (n in seq(1:14)){
	data<-read.table(paste('hmm.chr',n,'/',p,'_chr_',n,'_3rd_state.hmm.txt',sep=''),sep='\t',header=T)
	colnames(gap)<-colnames(data)
	mat<-rbind(mat,data,gap)	
}
n<-15
data<-read.table(paste('hmm.chr',n,'/',p,'_chr_',n,'_3rd_state.hmm.txt',sep=''),sep='\t',header=T)
mat<-rbind(mat,data)

p1<-unlist(strsplit(p,".",fixed=T))[1]
p1<-unlist(strsplit(p1,"pop",fixed=T))[2]
p2<-unlist(strsplit(p,".",fixed=T))[2]
p2<-unlist(strsplit(p2,"pop",fixed=T))[2]
mycolors<-c('orange',"darkgrey",'red')
mypch<-c(1,1,19)
mycex<-c(3,3,4)
if (args[1]=='Dsyl'){
plot(1,type="n",xlim=c(1,length(mat[,2])),ylim=c(0,0.8),xlab="",ylab=paste('Pop',p1,'.Pop',p2,sep=''),cex.lab=6,axes=F,mgp=c(10,4,0))
}
if (args[1]=='Dcar'){
plot(1,type="n",xlim=c(1,length(mat[,2])),ylim=c(0,1),xlab="",ylab=paste('Pop',p1,'.Pop',p2,sep=''),cex.lab=6,axes=F,mgp=c(10,4,0))
}


points(1:length(mat[,2]),mat[,5],col=mycolors[mat[,6]],pch=mypch[mat[,6]],cex=mycex[mat[,6]])
axis(1,labels=F,tck=0,lwd=6)
axis(2,labels=T,cex.axis=5,pos=-2000,mgp=c(3,2,0),lwd=6)

x<-paste('efst.',p,sep='')
mat2<-subset(mat,get(x)>0.005)
hist(mat2[,5],breaks=1000,main='',xlab="",ylab='',xlim=c(0,0.6),axes=F,col='black')
axis(1,labels=F,lwd=6)
axis(2,labels=T,cex.axis=5,lwd=6)

}

# ------------------------------------------------------------------------------------------------------
# Here I plot the other 6 (-1) PW within altitude
# ------------------------------------------------------------------------------------------------------

for (p in pops2[-length(pops2)]){
print(p)

mat<-c()
for (n in seq(1:14)){
	data<-read.table(paste('hmm.chr',n,'/',p,'_chr_',n,'_3rd_state.hmm.txt',sep=''),sep='\t',header=T)
	colnames(gap)<-colnames(data)
	mat<-rbind(mat,data,gap)	
}
n<-15
data<-read.table(paste('hmm.chr',n,'/',p,'_chr_',n,'_3rd_state.hmm.txt',sep=''),sep='\t',header=T)
mat<-rbind(mat,data)

p1<-unlist(strsplit(p,".",fixed=T))[1]
p1<-unlist(strsplit(p1,"pop",fixed=T))[2]
p2<-unlist(strsplit(p,".",fixed=T))[2]
p2<-unlist(strsplit(p2,"pop",fixed=T))[2]
mycolors<-c('orange',"darkgrey",'red')
mypch<-c(1,1,19)
mycex<-c(3,3,4)
if (args[1]=='Dsyl'){
plot(1,type="n",xlim=c(1,length(mat[,2])),ylim=c(0,0.8),xlab="",ylab=paste('Pop',p1,'.Pop',p2,sep=''),cex.lab=6,axes=F,mgp=c(10,4,0))
}
if (args[1]=='Dcar'){
plot(1,type="n",xlim=c(1,length(mat[,2])),ylim=c(0,1),xlab="",ylab=paste('Pop',p1,'.Pop',p2,sep=''),cex.lab=6,axes=F,mgp=c(10,4,0))
}
points(1:length(mat[,2]),mat[,5],col=mycolors[mat[,6]],pch=mypch[mat[,6]],cex=mycex[mat[,6]])
axis(1,labels=F,tck=0,lwd=6)
axis(2,labels=T,cex.axis=5,pos=-2000,mgp=c(3,2,0),lwd=6)

x<-paste('efst.',p,sep='')
mat2<-subset(mat,get(x)>0.005)
hist(mat2[,5],breaks=1000,main='',xlab="",ylab='',xlim=c(0,0.6),axes=F,col='black')
axis(1,labels=F,lwd=6)
axis(2,labels=T,cex.axis=5,lwd=6)

}

# ------------------------------------------------------------------------------------------------------
# Here I plot the last PW within altitude, with plot labels
# ------------------------------------------------------------------------------------------------------

p<-pops2[length(pops2)]

mat<-c()
for (n in seq(1:14)){
	data<-read.table(paste('hmm.chr',n,'/',p,'_chr_',n,'_3rd_state.hmm.txt',sep=''),sep='\t',header=T)
	colnames(gap)<-colnames(data)
	mat<-rbind(mat,data,gap)	
}
n<-15
data<-read.table(paste('hmm.chr',n,'/',p,'_chr_',n,'_3rd_state.hmm.txt',sep=''),sep='\t',header=T)
mat<-rbind(mat,data)

p1<-unlist(strsplit(p,".",fixed=T))[1]
p1<-unlist(strsplit(p1,"pop",fixed=T))[2]
p2<-unlist(strsplit(p,".",fixed=T))[2]
p2<-unlist(strsplit(p2,"pop",fixed=T))[2]
mycolors<-c('orange',"darkgrey",'red')
mypch<-c(1,1,19)
mycex<-c(3,3,4)
if (args[1]=='Dsyl'){
plot(1,type="n",xlim=c(1,length(mat[,2])),ylim=c(0,0.8),xlab="Chromosome",ylab=paste('Pop',p1,'.Pop',p2,sep=''),cex.lab=6,axes=F,mgp=c(10,4,0))
}
if (args[1]=='Dcar'){
plot(1,type="n",xlim=c(1,length(mat[,2])),ylim=c(0,1),xlab="Chromosome",ylab=paste('Pop',p1,'.Pop',p2,sep=''),cex.lab=6,axes=F,mgp=c(10,4,0))
}
points(1:length(mat[,2]),mat[,5],col=mycolors[mat[,6]],pch=mypch[mat[,6]],cex=mycex[mat[,6]])

if (args[1]=='Dsyl'){
xlab<-c(0,6000,16000,26000,37000,47000,58000,67000,76000,84500,94000,102500,111000,121000,134000,145000,170000)
}
if (args[1]=='Dcar'){
xlab<-c(0,5000,16000,27000,37000,47000,55000,63000,74000,87500,100000,112000,123000,134000,145000,155000,170000)
}
axis(1,labels=c('',1:15,''),at=xlab,lwd=6,tck=0,cex.axis=5,mgp=c(3,4,0))
axis(2,labels=T,cex.axis=5,pos=-2000,lwd=5,mgp=c(3,4,0))



x<-paste('efst.',p,sep='')
mat2<-subset(mat,get(x)>0.005)
hist(mat2[,5],breaks=1000,main='',xlab="Fst",ylab='',xlim=c(0,0.6),axes=F,col='black',mgp=c(10,4,0),cex.lab=6)
axis(1,labels=T,cex.axis=5,lwd=6,mgp=c(3,4,0))
axis(2,labels=T,cex.axis=5,lwd=6)



dev.off()




# ====================================================================================================
# Now plot the unmapped scaffolds
# ====================================================================================================


png(file=paste("combined_plots_unmapped_scaffs_Fst",".png",sep=""),width=5000,height=6000)
layout(matrix(c(1:30),15,2,byrow = TRUE), widths=c(10,1.5),heights=c(1,1),TRUE) 
par(mar=c(12,15,1,1))
par(oma=c(20,10,10,10))

# ------------------------------------------------------------------------------------------------------
#  1rst PW across altitude to write candidate numbers on top
# ------------------------------------------------------------------------------------------------------


for (p in pops1[1]){
print(p)

n<-16
data<-read.table(paste('hmm.chr',n,'/',p,'_chr_',n,'_3rd_state.hmm.txt',sep=''),sep='\t',header=T)

mat<-data
p1<-unlist(strsplit(p,".",fixed=T))[1]
p1<-unlist(strsplit(p1,"pop",fixed=T))[2]
p2<-unlist(strsplit(p,".",fixed=T))[2]
p2<-unlist(strsplit(p2,"pop",fixed=T))[2]
mycolors<-c('orange',"darkgrey",'red')
mypch<-c(1,1,19)
mycex<-c(3,3,4)
if (args[1]=='Dsyl'){
plot(1,type="n",xlim=c(1,length(mat[,2])),ylim=c(0,0.8),xlab="",ylab=paste('Pop',p1,'.Pop',p2,sep=''),cex.lab=6,axes=F,mgp=c(10,4,0))
mtext(c(14:33),at=c(1990,7289,23297,30133,42626,67000,70514,83080,86320,96327,99488,104134,108426,126168,133401,162194,171322,176657,189222,192731),cex=3,side=3,col='red')
}
if (args[1]=='Dcar'){
plot(1,type="n",xlim=c(1,length(mat[,2])),ylim=c(0,1),xlab="",ylab=paste('Pop',p1,'.Pop',p2,sep=''),cex.lab=6,axes=F,mgp=c(10,4,0))
mtext(c(7:27),at=c(35065,47380,49305,54500,57000,59698,62444,82000,84500,87000,89500,92000,95075,101229,105296,108500,111500,126635,143182,150500,161127),cex=3,side=3,col='red')
}
points(1:length(mat[,2]),mat[,5],col=mycolors[mat[,6]],pch=mypch[mat[,6]],cex=mycex[mat[,6]])
axis(1,labels=F,tck=0,lwd=6)
axis(2,labels=T,cex.axis=5,pos=-2000,mgp=c(3,2,0),lwd=6)



x<-paste('efst.',p,sep='')
mat2<-subset(mat,get(x)>0.005)
hist(mat2[,5],breaks=1000,main='',xlab="",ylab='',xlim=c(0,0.6),axes=F,col='black')
axis(1,labels=F,lwd=6)
axis(2,labels=T,cex.axis=5,lwd=6)

}

# ------------------------------------------------------------------------------------------------------
#  rest of the PW 
# ------------------------------------------------------------------------------------------------------

for (p in c(pops1[-1],pops2)){
print(p)

n<-16
data<-read.table(paste('hmm.chr',n,'/',p,'_chr_',n,'_3rd_state.hmm.txt',sep=''),sep='\t',header=T)

mat<-data
p1<-unlist(strsplit(p,".",fixed=T))[1]
p1<-unlist(strsplit(p1,"pop",fixed=T))[2]
p2<-unlist(strsplit(p,".",fixed=T))[2]
p2<-unlist(strsplit(p2,"pop",fixed=T))[2]
mycolors<-c('orange',"darkgrey",'red')
mypch<-c(1,1,19)
mycex<-c(3,3,4)
if (args[1]=='Dsyl'){
plot(1,type="n",xlim=c(1,length(mat[,2])),ylim=c(0,0.8),xlab="",ylab=paste('Pop',p1,'.Pop',p2,sep=''),cex.lab=6,axes=F,mgp=c(10,4,0))
}
if (args[1]=='Dcar'){
plot(1,type="n",xlim=c(1,length(mat[,2])),ylim=c(0,1),xlab="",ylab=paste('Pop',p1,'.Pop',p2,sep=''),cex.lab=6,axes=F,mgp=c(10,4,0))
}
points(1:length(mat[,2]),mat[,5],col=mycolors[mat[,6]],pch=mypch[mat[,6]],cex=mycex[mat[,6]])
axis(1,labels=F,tck=0,lwd=6)
axis(2,labels=T,cex.axis=5,pos=-2000,mgp=c(3,2,0),lwd=6)

x<-paste('efst.',p,sep='')
mat2<-subset(mat,get(x)>0.005)
hist(mat2[,5],breaks=1000,main='',xlab="",ylab='',xlim=c(0,0.6),axes=F,col='black')
axis(1,labels=F,lwd=6)
axis(2,labels=T,cex.axis=5,lwd=6)

}

dev.off()

}




# ====================================================================================================
# This is to make the PDF instead of single pops
# ====================================================================================================

 if (args[4] == 'pdf'){
 for (p in pops1){
 	 	
pdf(paste(p,"_plots_Fst",".pdf",sep=""),onefile=T,width = 60, height = 15)
layout(matrix(c(1:2),1,2,byrow = TRUE), widths=c(10,1.5),heights=c(2,2),TRUE) 
par(mar=c(10,10,1,1))
par(oma=c(20,10,10,10))

# ------------------------------------------------------------------------------------------------------
# Here I plot the chrs
# ------------------------------------------------------------------------------------------------------

mat<-c()
for (n in seq(1:14)){
	data<-read.table(paste('hmm.chr',n,'/',p,'_chr_',n,'_3rd_state.hmm.txt',sep=''),sep='\t',header=T)
	colnames(gap)<-colnames(data)
	mat<-rbind(mat,data,gap)	
}
n<-15
data<-read.table(paste('hmm.chr',n,'/',p,'_chr_',n,'_3rd_state.hmm.txt',sep=''),sep='\t',header=T)
mat<-rbind(mat,data)

p1<-unlist(strsplit(p,".",fixed=T))[1]
p1<-unlist(strsplit(p1,"pop",fixed=T))[2]
p2<-unlist(strsplit(p,".",fixed=T))[2]
p2<-unlist(strsplit(p2,"pop",fixed=T))[2]
mycolors<-c('orange',"darkgrey",'red')
mypch<-c(1,1,19)
mycex<-c(1,1,2)
if (args[1]=='Dsyl'){
plot(1,type="n",xlim=c(1,length(mat[,2])),ylim=c(0,0.8),xlab="Chromosome",ylab=paste('Pop',p1,'.Pop',p2,sep=''),cex.lab=4,axes=F,mgp=c(7,3,0))
}
if (args[1]=='Dcar'){
plot(1,type="n",xlim=c(1,length(mat[,2])),ylim=c(0,1),xlab="Chromosome",ylab=paste('Pop',p1,'.Pop',p2,sep=''),cex.lab=4,axes=F,mgp=c(7,3,0))
}
points(1:length(mat[,2]),mat[,5],col=mycolors[mat[,6]],pch=mypch[mat[,6]],cex=mycex[mat[,6]])



if (args[1]=='Dsyl'){
xlab<-c(0,6000,16000,26000,37000,47000,58000,67000,76000,84500,94000,102500,111000,121000,134000,145000,170000)
mtext(c(1:13),at=c(6500,9500,30500,36000,38000,39200,75500,94000,117500,123500,131500,136500,142000),cex=3,side=3,col='red')
}
if (args[1]=='Dcar'){
xlab<-c(0,5000,16000,27000,37000,47000,55000,63000,74000,87500,100000,112000,123000,134000,145000,155000,170000)
mtext(c(1:6),at=c(67000,73000,75500,77000,90500,153000),cex=3,side=3,col='red')
}
axis(1,labels=c('',1:15,''),at=xlab,lwd=5,tck=0,cex.axis=3,mgp=c(7,3,0))
axis(2,labels=T,cex.axis=3,pos=-2000,lwd=5,mgp=c(3,2,0))


x<-paste('efst.',p,sep='')
mat2<-subset(mat,get(x)>0.005)
hist(mat2[,5],breaks=1000,main='',xlab="Fst",ylab='',xlim=c(0,0.6),axes=F,col='black',cex.lab=4,mgp=c(7,3,0))
axis(1,labels=T,cex.axis=3,lwd=6,mgp=c(7,3,0))
axis(2,labels=T,cex.axis=3,lwd=6)

}

# ------------------------------------------------------------------------------------------------------
# Here I plot the unmapped scaffolds
# ------------------------------------------------------------------------------------------------------

n<-16
data<-read.table(paste('hmm.chr',n,'/',p,'_chr_',n,'_3rd_state.hmm.txt',sep=''),sep='\t',header=T)
mat<-data
mat<-data.frame(mat)
print(head(mat))



### ----------------------------------------
#   This is to find the positions of the candidates  
x<-paste('state.',p,sep='')
cand<-subset(mat,get(x)=='3')
cand<-data.frame(cand)
cand.scaffs<-unique(cand[,1])
for (scaff in cand.scaffs){
	print(subset(cand,chr==scaff)[1,][1])
}
##-------------------------------------------

p1<-unlist(strsplit(p,".",fixed=T))[1]
p1<-unlist(strsplit(p1,"pop",fixed=T))[2]
p2<-unlist(strsplit(p,".",fixed=T))[2]
p2<-unlist(strsplit(p2,"pop",fixed=T))[2]
mycolors<-c('orange',"darkgrey",'red')
mypch<-c(1,1,19)
mycex<-c(1,1,2)

if (args[1]=='Dsyl'){
plot(1,type="n",xlim=c(1,length(mat[,2])),ylim=c(0,0.8),xlab="",ylab=paste('Pop',p1,'.Pop',p2,sep=''),cex.lab=4,axes=F,mgp=c(7,3,0))
}
if (args[1]=='Dcar'){
plot(1,type="n",xlim=c(1,length(mat[,2])),ylim=c(0,1),xlab="",ylab=paste('Pop',p1,'.Pop',p2,sep=''),cex.lab=4,axes=F,mgp=c(7,3,0))
}
points(1:length(mat[,2]),mat[,5],col=mycolors[mat[,6]],pch=mypch[mat[,6]],cex=mycex[mat[,6]])

if (args[1]=='Dsyl'){
mtext(c(14:33),at=c(1990,7289,23297,30133,42626,67000,70514,83080,86320,96327,99488,104134,108426,126168,133401,162194,171322,176657,189222,192731),cex=3,side=3,col='red')
}
if (args[1]=='Dcar'){
mtext(c(7:27),at=c(35065,47380,49305,54500,57000,59698,62444,82000,84500,87000,89500,92000,95075,101229,105296,108500,111500,126635,143182,150500,161127),cex=3,side=3,col='red')
}

axis(2,labels=T,cex.axis=3,pos=-2000,lwd=5,mgp=c(3,2,0))


x<-paste('efst.',p,sep='')
mat2<-subset(mat,get(x)>0.005)
hist(mat2[,5],breaks=1000,main='',xlab="Fst",ylab='',xlim=c(0,0.6),axes=F,col='black',cex.lab=4,mgp=c(7,3,0))
axis(1,labels=T,cex.axis=3,lwd=6,mgp=c(7,3,0))
axis(2,labels=T,cex.axis=3,lwd=6)


dev.off()

}







