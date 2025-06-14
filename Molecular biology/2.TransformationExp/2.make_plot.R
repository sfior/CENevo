
traits<-c('TLN_SD','days_till_bolting_SD','TLN_LD','days_till_bolting_LD')


mat_mean<-c()
mat_sd<-c()
mat_se<-c()
mat_ci<-c()

for (t in traits){
print(t)
data<-read.table(paste(t,'.stats.txt',sep=''),header=T,sep='\t',na.string='NA',fill=T)

mat_mean<-cbind(mat_mean,data$x)
mat_sd<-cbind(mat_sd,data$sd)
mat_se<-cbind(mat_se,data$se)
mat_ci<-cbind(mat_ci,data$ci)
}

rownames(mat_mean)<-c('Col','High','Low')
rownames(mat_sd)<-c('Col','High','Low')
rownames(mat_se)<-c('Col','High','Low')
rownames(mat_ci)<-c('Col','High','Low')

colnames(mat_mean)<-c('TLN_SD','days_till_bolting_SD','TLN_LD','days_till_bolting_LD')
colnames(mat_sd)<-c('TLN_SD','days_till_bolting_SD','TLN_LD','days_till_bolting_LD')
colnames(mat_se)<-c('TLN_SD','days_till_bolting_SD','TLN_LD','days_till_bolting_LD')
colnames(mat_ci)<-c('TLN_SD','days_till_bolting_SD','TLN_LD','days_till_bolting_LD')

error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
}

trait_barplot = barplot(mat_mean , beside=T , legend.text=T,col=c("green" , "blue", "red"), ylim=c(0,80) , ylab="yyy",las=2)
error.bar(trait_barplot,mat_mean, mat_se)


