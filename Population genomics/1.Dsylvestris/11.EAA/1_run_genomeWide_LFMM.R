##### LFMM with D. Sylvestris ######

### load libraries
library(matrixStats)
library(lfmm)
### Set directory
setwd("./")

### Read environmental data
env <- read.table("env_pops.txt", header=T, sep="\t")
env <- env[,4:ncol(env)]
env <- as.data.frame(env)

### Set the variables
# Env variable column in table
Ecol<-1 # Elevation

# Env variable 
Evar<-colnames(env)[Ecol]
print(Evar)

### Read allele frequency data
mafs <- read.table("AllPops_innerMerged.AF.AFcolumnsOnly.gz", header=T)
#mafs<-mafs[1:1000,]   # subset
mafp <- mafs[,-1]
rownames(mafp) <- mafs$locus

### Filter allele frequencies
#mafmax <- rowMaxs(as.matrix(mafp))
#mafmin <- rowMins(as.matrix(mafp))
#mafp <- mafp[-c(which(mafmax<0.075),which(mafmin>0.925)),]
mafs <- mafp
remove(mafp)
#head(mafs)

### LFMM
Y1 <- t(mafs)
K <- 2
gif.lfmm <- matrix(nrow = 1, ncol = K)
rownames(gif.lfmm) <- colnames(Evar)
colnames(gif.lfmm) <- paste("K", 1:K, sep = "")
fdr.thres <- c(0.05)

for(i in 1:K){
  
  mod.lfmm <- lfmm_ridge(Y=Y1, X=env[,Ecol], K=i)
  stats.lfmm <- lfmm_test(Y=Y1, X=env[,Ecol], lfmm=mod.lfmm, calibrate="gif")
  
  results.lfmm <- matrix(nrow=NCOL(Y1), ncol=5)
  rownames(results.lfmm) <- colnames(Y1)
  colnames(results.lfmm) <- c("snpid","zscore","pvalue","adjpvalue","beta")
  
  gif.lfmm[1,i] <- round(stats.lfmm$gif,2)
  results.lfmm[,"snpid"] <- colnames(Y1)
  results.lfmm[,"zscore"] <- stats.lfmm$score
  results.lfmm[,"pvalue"] <- stats.lfmm$pvalue
  results.lfmm[,"adjpvalue"] <- stats.lfmm$calibrated.pvalue
  results.lfmm[,"beta"] <- stats.lfmm$B
  
  
  # Write all results
  write.table(results.lfmm, paste("LFMM_results_",Evar,"_K", i, ".txt", sep=""), row.names=F, col.names=T, quote = F, sep='\t') #All results
  
  # FDR correction
  results.lfmm <- as.data.frame(results.lfmm)
  results.lfmm.sign <- results.lfmm[which( p.adjust(as.vector(as.numeric(results.lfmm[,"adjpvalue"])), method="fdr", n=NROW(results.lfmm)) <= fdr.thres),]
#  if (nrow(as.data.frame(results.lfmm.sign))>2){
    write.table(results.lfmm.sign, paste("LFMM_results_",Evar,"_K", i, "q", fdr.thres, ".txt", sep=""), sep="\t", row.names=F, col.names=T, quote=F)
    # sort according to adjpvalue
    results.lfmm.sign <- results.lfmm.sign[order(as.numeric(results.lfmm.sign[,"adjpvalue"]),decreasing=FALSE),]
    write.table(results.lfmm.sign, paste("LFMM_results_",Evar,"_K", i, "q", fdr.thres, "_sorted.txt", sep=""), sep="\t", row.names=F, col.names=T, quote=F)
#  }
}
write.table(gif.lfmm, paste("GIF","_K1-", K,".txt", sep=""), sep=";", row.names=F, col.names=T)



