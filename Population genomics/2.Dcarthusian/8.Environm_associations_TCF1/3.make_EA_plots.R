
# This script makes plots of p-values for the associations between allele frequencies and an environmental variable (partial Mantel test)
# It outputs PNG plots for each variable

# This script is run from the command line with Rscript.
# Requires:
# 1: scaffold names used in input files
# 2: list of env. variable names used in input files

# Usage:
# Rscript make_EA_plots.R scaffold65 sfroyy_100,srad_average 


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

autoloess <- function(fit, span=c(as.numeric(0.001), as.numeric(0.01))) {
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




args <- commandArgs(trailingOnly = TRUE)
if(length(args)==6){
    print("Dataset will not be subset")
}

scaffs<-strsplit(args[1],",")[[1]]
vars<-strsplit(args[2],",")[[1]]  


for (s in scaffs){
	
	for (v in vars){
		data<-read.table(paste(s,'_',v,'.txt',sep=''),header=T,sep='\t')
		b<-data[,2]
		p1<-log10(data[,9])
		pv<-cbind(b,p1)
		set.seed(62485)
		basic <- loess(p1 ~ b)
		good <- autoloess(basic)
		opt_span<-good$pars$span
		print(opt_span)
		png(paste(s,'_',v,".png",sep=""),width = 4, height = 3, units = 'in', res = 300)
		plot(pv,cex=0.01,col="lightgrey",xlab='position',ylab="log(pvalue)",cex.axis=0.5,cex.lab=0.7,tck=-0.02,mgp=c(1,0.1,0))
		#title(ylab="log(pvalue)", mgp=c(1,0.3,0), cex.lab=0.7 )
		#legend("bottomleft",paste("span=",round(opt_span,digits=3),sep=""),lwd=1,cex=0.5,bty='n')
		lim <- par("usr")
		peakx1<-156250
		peakx2<-166250
		rect(peakx1, lim[3]-1, peakx2, lim[4]+1, border = "khaki1", col = "khaki1", density= 30)

		lines(b,good$fitted,col=1,lwd=1)
		dev.off()		
	}	
}