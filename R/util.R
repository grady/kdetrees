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


##' Reorder objects in a kdetrees object based on current density scores.
##'
##' @param x kdetrees object
##' @param decreasing order of sorting
##' @param ... additional parameters passed to \code{\link{order}}
##' @return a new kdetrees object
##' @export
##' @author Grady Weyenberg
sort.kdetrees <- function(x, decreasing = FALSE, ...){
  oi <- order(x$density, decreasing = decreasing,...) #ordered indices

  x$km <- x$km[oi,oi]
  x$density <- x$density[oi]
  x$i <- seq_along(x$i)
  x$trees <- x$trees[oi]
  x$outliers <- x$trees[x$i]
  x
}

write.kdetrees <- function(x,dir=".",trees=TRUE,...){
  if(!file.exists(dir)) dir.create(dir,recursive=TRUE)
  oldDir <- setwd(dir)
  on.exit(setwd(oldDir))
  ggsave("scatter.pdf", plot(x))
  ggsave("histogram.pdf",hist(x))
  df <- data.frame(score=x$density,
                   outlier=(seq_along(x$density) %in% x$i),
                   file=names(x$trees))
  if(trees) df$tree=write.tree(x$trees)
  write.table(df, "scores.txt", row.names=FALSE)
  if(length(x$outliers) > 0){
      write.tree(x$outliers, "outliers.tre", tree.names=TRUE)
      write.tree(x$trees[-x$i], "non-outlier.tre", tree.names=TRUE)
  } else {
      write.tree(x$trees, "non-outlier.tre", tree.names=TRUE)
  }
  invisible(x)
}



## Suppress the NOTES from R CMD check about undefined variables (ggplot calls)
globalVariables(c("outlier","index"))


zero.leaf.edges <- function(x){
    x <- unroot(x)
    n <- Ntip(x)
    is.leaf <- function(z) {any(z %in% 1:n)}
    i <- apply(x$edge,1, is.leaf)
    x$edge.length[i] <- 0.0
    x
}

##' Read Newick trees from a set of files.
##'
##' Read a set of files containing individual trees and return a
##' multiPhylo object with all trees in the files. The set of files to
##' be read can either be passed explicitly as a character vector, or
##' through a combination of a set of directories and a file extension
##' which is used to construct a wildcard expression which will match
##' files with the extension you specify..
##'
##' The behavior is determined by the extension parameter: If it is
##' NULL, then the paths supplied are passed unaltered to the
##' ape::read.tree function. If the parameter is not null, then the
##' paths are altered by calling Sys.glob with wildcard expression
##' "paths*extension".
##' 
##' @param paths If estension is NULL, a list of files containing
##' trees. Otherwise, the prefix for a wildcard expression.
##' @param extension If NULL, the paths are passed unaltered to
##' read.tree. Otherwise, a suffix (file extension) for the wildcard
##' expression.
##' @param use.file.names Set the tree names using the file
##' names. Will only work properly if each file contains a single
##' tree.
##' @return A multiPhylo object.
##' @author Grady Weyenberg
##' @export
load.trees <- function(paths,extension=".tre",use.file.names=TRUE){
    if(!is.null(extension))
        paths <- Sys.glob(paste(paths,"*",extension,sep=""))
    trees <- lapply(paths, read.tree)
    trees <- lapply(trees, function(x){x$node.label <- NULL;x})
    trees <- do.call(c,trees)
    if(use.file.names){
        if(length(paths) == length(trees))
            names(trees) <- basename(paths)
        else
            warning("Found ",length(paths)," files and ",length(trees),
                    "trees. Ignoring use.file.names.")
    }
    trees
}

plot.multiPhylo <- function (x, layout = 1, ...) {
    layout(matrix(1:layout, ceiling(sqrt(layout)), byrow = TRUE))
    ## if (!devAskNewPage() && interactive()) {
    ##     devAskNewPage(TRUE)
    ##     on.exit(devAskNewPage(FALSE))
    ## }
    for (i in 1:length(x)){
        plot(x[[i]], main=names(x)[i],...)
    }
}
