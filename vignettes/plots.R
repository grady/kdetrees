load("sim/res.Rdata")
load("sim/res2.Rdata")

library(ggplot2)
library(plyr)
library(BatchExperiments)

res2 <- ddply(res, getResultVars(res),summarise,hit=mean(hit),type1=mean(type1))
levels(res2$distance) <- c("kde-dissimilarity","kde-geodesic","pmc-topological","pmc-edges")
levels(res2$sp.file) <- c("single","2-mixture","5-mixture")

ggplot(subset(res2,!(distance=="pmc-topological")), aes(x=type1,y=hit,col=factor(neff)))+geom_line()+facet_grid(distance~sp.file)+geom_abline(intercept=0,slope=1,lty=3,alpha=0.2)+xlim(0,1)+ylim(0,1)+labs(title="ROC Plots",x="FPR",y="TPR")+scale_color_grey(expression(N[e]),end=0.75,start=0)



ggplot(subset(res2,(k==1.25 & algo=="kdetrees")|(k==0.125 & algo=="pmcoa")), aes(x=neff,y=hit,col=algo,pch=distance))+geom_point()+facet_grid(sp.file~.)+geom_smooth(se=FALSE)+ggtitle("k selected based on Type1 rate") #+geom_smooth(aes(y=type1),lty=2)#+geom_abline(intercept=0.05,slope=0,lty=2)

## 2 > k > -1.3 for pmcoa
## 2 > k for kde not single up to 3 for single
## k > -1.5 for kde not sp5 up to -2 for sp5


ggplot(res, aes(x=type1,y=hit,col=neff,group=factor(neff)))+stat_summary(fun.data=mean_cl_normal,geom="errorbar")+facet_grid(distance~sp.file)+geom_abline(intercept=0,slope=1,alpha=0.2)+ylim(0,1)+

ggplot(res, aes(x=type1,y=hit,col=neff,group=factor(neff)))+stat_summary(fun.y=mean,geom="line")+facet_grid(distance~sp.file)+geom_abline(intercept=0,slope=1,alpha=0.2)+ylim(0,1)

ggplot(res2, aes(x=neff,y=hit,lty=distance,col=(k),group=factor(k)))+geom_point()+facet_grid(.~sp.file)+geom_smooth(se=FALSE)+scale_x_log10(breaks=5)

