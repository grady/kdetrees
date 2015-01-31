### Copyright (C) 2014 -- Grady Weyenberg ###
##' Summarize a kdetrees object in human-readable form.
##'
##' Pretty-prints the results of a kdetrees ananlysis to console.
##'
##' @param x object to be printed
##' @param ... unused, required for generic compatability
##' @return invisible(x)
##' @author Grady Weyenberg
##' @method print kdetrees
##' @export
print.kdetrees <- function(x,...){
  outliers <- if(is.null(names(x$outliers))) x$i else names(x$outliers)
  cat("Call: ")
  print(attr(x,"call"))
  cat("Density estimates:\n")
  print(summary(x$density))
  cat("Cutoff: ",attr(x,"c"),"\n")
  cat("\nOutliers detected:\n")
  print(outliers,quote=FALSE)
  invisible(x)
}

##' Convert kdetrees object to data.frame
##'
##' Converts a kdetrees object to a data.frame suitable for saving as
##' output. It contains the density estimates for each tree, a Boolean
##' value indicating if the tree was selected as an outlier, and
##' optionally the newick string corresponding to the tree.
##' 
##' @param x kdetrees object to be converted
##' @param row.names ignored
##' @param optional ignored
##' @param trees If given the original list of trees, will convert to
##' newick and add a column to the output
##' @param ... unused
##' @return a data.frame
##' @author Grady Weyenberg
##' @method as.data.frame kdetrees
##' @export
##' @examples
##' result <- kdetrees(apicomplexa)
##' as.data.frame(result)
as.data.frame.kdetrees <- function(x, row.names, optional, trees=NULL, ...){
    out <- data.frame(density = x$density,
                      outlier=names(x$density) %in% names(x$outliers),
                      stringsAsFactors=FALSE)
    if (inherits(trees,"multiPhylo")) out$newick <- write.tree(trees,digits=4L)
    out
}
  
##' Plot the unnormalized density estimates for each tree.
##'
##' @param x kdetrees object to be plotted
##' @param ... additional arguments passed to ggplot
##' @return a ggplot object
##' @author Grady Weyenberg
##' @export
##' @method plot kdetrees
##' @examples
##' result <- kdetrees(apicomplexa)
##' plot(result)
plot.kdetrees <- function(x,...){
  df <- with(x,data.frame(density=unname(density),
                          index=seq_along(density),
                          outlier=seq_along(density) %in% i))
  ylab <- "Non-normalized Density"
  xlab <- "Tree Index"
  main <- paste(length(x$outliers),"Outliers Removed")
  
  ggplot(df,aes(index,density,color=outlier),...) + geom_point() +
    labs(title=main, x=xlab,y=ylab) + theme(legend.position="top")
}


##' Create a histogram of tree density estimates
##' 
##' @param x kdetrees object to plot
##' @param ... additional arguments passed to ggplot
##' @return a ggplot object
##' @author Grady Weyenberg
##' @export
##' @method hist kdetrees
##' @examples
##' result <- kdetrees(apicomplexa)
##' hist(result)
hist.kdetrees <- function(x,...){
  df <- with(x,data.frame(density=unname(density),
                          index=seq_along(density),
                          outlier=seq_along(density) %in% i))
  bw <- with(x,diff(range(density))/nclass.FD(density))
  main <- paste("Histogram of Estimates:",length(x$outliers),"Outliers Removed")
  xlab <- "Non-normalized Density"
  ylab <- "Count"

  ggplot(df,aes(density,fill=outlier),...) + geom_histogram(binwidth=bw) +
    labs(title=main,x=xlab,y=ylab) + theme(legend.position="top")
}

##' Extract vector of internal branch lengths
##'
##' Take a phylo object and extract a vector containing the internal
##' edge lengths, in no particular order.
##' @param x a phylo object
##' @return a vector with the internal edge lengths
##' @author Grady Weyenberg
internal.blens <- function(x){
    x <- unroot(x)
    n <- Ntip(x)
    f <- function(z) {any(z %in% 1:n)}
    i <- apply(x$edge,1,f)
    sort(x$edge.length[!i])
}

##' count semi-labeled rooted trees on n taxa
##'
##' counts the number of topologies for rooted trees on n taxa
##' (NOTE) Candidate for optimization.
##' 
##' @param n number of taxa, must be at least 3
##' @return integer number of topologies for a rooted phylogeny
##' @author Grady Weyenberg
Northant <- function(n){
    n <- trunc(n)
    if (any(n<3)) stop("n must be greater than 2")
    x <- 2L * n - 3L
    res <- rep(1L,length(n))
    while (any(x > 1L)) {
        res <- res * x
        x <- x - 2L
        x[x<1L] <- 1L
    }
    res
}

find.const <- function(blens,bw,nOrth){
    blens <- sort(blens)
    d <- length(blens)
    bw <- diag(d) * bw
    c0 <- 1/sqrt((2 * pi)^d * det(bw)) #gaussian const at center
    p1 <- hgm.ncorthant(bw, blens)
    p0 <- hgm.ncorthant(bw, -blens)

    blens[1:2] <- -blens[1:2]
    p2 <- hgm.ncorthant(bw, blens)

    ll <- (p1 + p0*(nOrth-1))/c0
    ul <- (p1 + p2*(nOrth-1))/c0
    c(ll=ll,ul=ul)
}
## need a function that returns for each
bhv.consts <- function(trees,bw){
    x <- sapply(trees,internal.blens)
    if(!is.matrix(x)) stop ("different Ntips in trees")
    n <- Northant(Ntip(trees[[1]])) #number of orthants
    d <- nrow(x) #dim of orthants
    N <- ncol(x) #number of trees
    res <- numeric(N) #result obj

    if ( length(bw) == 1L ) bw <- rep(bw,N)

    for (i in 1:N){
        bm <- diag(d) * bw[i]
        c0 <- 1/sqrt((2 * pi)^d * det(bm)) #gaussian const at center
        p1 <- hgm.ncorthant(bm, x[,i])
        p0 <- hgm.ncorthant(bm, -x[,1])
        y <- x[,i]
        y[1] <- -y[1]
        p2 <- hgm.ncorthant(bm,y)
        
        res[i] <- (p1 + p0*(n-1))/c0
#        res[i] <- (p1 + p2 * 3 + p0 * (n-4))/c0
    }
#    print(x)
    res
}
## Suppress the NOTES from R CMD check about undefined variables (ggplot calls)
globalVariables(c("outlier","index"))

