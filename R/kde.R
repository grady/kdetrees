##' Identify discordant trees in a sample
##'
##' If bw is a single number, it will be used as a single constant
##' bandwidth. It can also be a vector, in which case it will be used
##' as variable bandwidths. Finally, if it is a list, the list will be
##' passed as arguments to the bw.nn adaptive bandwith function.
##' 
##' @param trees multiPhylo object
##' @param n number of outliers to report
##' @param bw see Details
##' @param ... additional arguments for dist.diss
##' @return a kdetrees object; list(density,outliers,bandwidth)
##' @author Grady Weyenberg
##' @export
kdetrees <- function(trees,n=ceiling(0.05*length(trees)),bw=list(),...){
  dm <- dist.diss(trees,...)
  if(is.list(bw)) bw <- do.call(bw.nn,c(list(dm),bw))
  km <- normkern(dm,bw)
  i <- which.min(estimate(km))
  while(length(i) < n){
    j <- which.min(estimate(km[-i,-i]))
    j[1] <- match(names(j),rownames(km))
    i <- c(i,j)
  }
  est <- estimate(km,i)
  out <- list(density=est,outliers=i,bandwidth=bw)
  class(out) <- "kdetrees"
  out
}

##' Plot the unnormalized density estimates for each tree.
##'
##' If ggplot2 graphics are used, ... is a sequence of additional
##' ggplot2 directives to be added to the end of the ggplot
##' construction. If base graphics are used, ... are additional
##' arguments passed to plot command.
##' @param x kdetrees object to be plotted
##' @param ... see details
##' @param ggplot Logical: ggplot2 or base graphics
##' @return either a ggplot object or NULL
##' @author Grady Weyenberg
##' @export
##' @method plot kdetrees
##' @examples
##' fit <- kdetrees(apicomplexa,12,use.blen=TRUE)
##' plot(fit)
plot.kdetrees <- function(x,...,ggplot=("ggplot2" %in% installed.packages())){
  df <- with(x,data.frame(density=unname(density),
                          tree=names(density),
                          index=seq_along(density),
                          outlier=seq_along(density)%in%outliers))
  ylab <- "Non-normalized Density"
  xlab <- "Tree Index"
  main <- paste(length(x$outliers),"Outliers Removed")
  
  if(ggplot)
    ggreduce(ggplot(df,aes(index,density,color=outlier)),
             geom_point(), labs(title=main, x=xlab,y=ylab),
             theme(legend.position="top"),...)
  else
    plot(density~index,data=df,pch=as.numeric(outlier)+1L,ylab=ylab,xlab=xlab,main=main)
}


##' Create a histogram of tree density estimates
##' 
##' If ggplot2 graphics are used, ... is a sequence of additional
##' ggplot2 directives to be added to the end of the ggplot
##' construction. If base graphics are used, ... are additional
##' arguments passed to plot command.
##' @param x kdetrees object to plot
##' @param ... see details
##' @param ggplot Logical: ggplot2 or base graphics
##' @return either a ggplot or histogram object
##' @author Grady Weyenberg
##' @export
##' @method hist kdetrees
##' @examples
##' fit <- kdetrees(apicomplexa,12,use.blen=TRUE)
##' plot(fit)
hist.kdetrees <- function(x,...,ggplot=("ggplot2" %in% installed.packages())){
  df <- with(x,data.frame(density=unname(density),
                          tree=names(density),
                          index=seq_along(density),
                          outlier=seq_along(density)%in%outliers))
  bw <- with(x,diff(range(density))/nclass.FD(density))
  main <- paste("Histogram of Estimates:",length(x$outliers),"Outliers Removed")
  xlab <- "Non-normalized Density"
  ylab <- "Count"

  if(ggplot)
    ggreduce(ggplot(df,aes(density,fill=outlier)),
             geom_histogram(binwidth=bw), labs(title=main,x=xlab,y=ylab),
             theme(legend.position="top"),...)
  else
    hist(df$density,main=main,xlab=xlab,ylab=ylab)
}
  

##' Makes ggplot work more like usual R. Uses ',' instead of '+'.
##' @param ... sequence of ggplot2 commands 
ggreduce <- function(...) eval.parent(substitute(Reduce("+",list(...))))

##' estimate densities from kernel matrix
##'
##' @param x matrix of kernel contributions
##' @param i vector of columns to exclude from calculation
##' @return vector of density estimates for each tree
##' @author Grady Weyenberg
estimate <- function(x,i=integer()){
  if(length(i) > 0)
    rowSums(x[,-i]) - (!(1:nrow(x) %in% i)) * diag(x)
  else
    rowSums(x) - diag(x)
}


