##' Identify discordant trees in a sample
##'
##' If bw is a single number, it will be used as a single constant
##' bandwidth. It can also be a vector, in which case it will be used
##' as variable bandwidths. Finally, if it is a list, the list will be
##' passed as arguments to the bw.nn adaptive bandwith function.
##'
##' ... Is passed to either \code{distory::dist.multiPhylo} or
##' \code{dist.diss}, as appropriate. See the help for these functions
##' for more details.
##' 
##' @param trees multiPhylo object
##' @param k IQR multiplier for outlier detection
##' @param distance Select "geodesic" or "dissimilarity" distance
##' calculation method
##' @param bw see Details
##' @param greedy greedy outlier detection?
##' @param ... additional arguments for distance calculation function, see details
##' @return a kdetrees object; list(density,outliers)
##' @author Grady Weyenberg
##' @export
##' @examples
##' kdeobj <- kdetrees(apicomplexa)
##' print(kdeobj)
##' kdeobj$outliers
##'
##' kdetrees(apicomplexa, k=2.0, distance="dissimilarity",use.blen=TRUE)
kdetrees <- function(trees,k=1.5,distance=c("geodesic","dissimilarity"),bw=list(),greedy=FALSE,...){
  distance <- match.arg(distance)
  dm <- switch(distance, geodesic = as.matrix(dist.multiPhylo(trees,...)),
               dissimilarity = as.matrix(dist.diss(trees,...)))
  dimnames(dm) <- list(names(trees),names(trees))

  cutoff <- function(x, c = 1.5){
    qs <- quantile(x, c(0.25,0.75))
    unname(diff(qs) * -c + qs[1])
  }
  
  if(is.list(bw)) bw <- do.call(bw.nn,c(list(dm),bw))
  km <- normkern(dm,bw)
  x <- estimate(km)
  c <- cutoff(x, k)
  i <- which( x < c )
  if (greedy) {
    while(TRUE){
      i <- which( x < c )
      if (length(i) < 1) break
      x <- estimate(km,i)
      c2 <- cutoff(x[-i], k)
      if(is.na(c2)) browser()
      if(c2 > c) c <- c2 else break
    }
  }
  structure(list(density=x, i=i, outliers=trees[i]), class="kdetrees",
            call=match.call(), c=c)
}
##' do complete analysis in one call
##' @param file newick file with trees
##' @param outgroup if a character, reroot all trees with this species as outgroup
##' @param ... additional parameters for kdetrees
##' @param tree.file write outlier trees in newick format to this file
##' @param csv.file write density results to this file
##' @param plot.file print scatterplot of results to this file
##' @param hist.file print histogram of density estimates to this file
##' @return results of kdetrees call
##' @author Grady Weyenberg
##' @export
kdetrees.complete <- function(file, outgroup=NULL,...,tree.file="outliers.tre",
                              csv.file="results.csv",plot.file="plot.png",
                              hist.file="hist.png"){
  trees <- read.tree(file)
  if (is.null(names(trees))) names(trees) <- paste("tree",seq_along(trees),sep="")
  if (!inherits(trees,"multiPhylo")) stop("Could not read tree file")
  if (is.character(outgroup)) {
    trees <- lapply(trees,root,outgroup,resolve.root=TRUE)
    trees <- lapply(trees,"[<-","node.label", NULL)
    class(trees) <- "multiPhylo"
  }
  
  res <- kdetrees(trees,...)
  browser()
  if (is.character(plot.file)) ggsave(plot.file,plot(res))
  if (is.character(hist.file)) ggsave(hist.file,hist(res))
  if (is.character(csv.file)) write.csv(as.data.frame(res),csv.file)
  if (is.character(tree.file) && length(res$outliers) > 0)
    write.tree(res$outliers, tree.file, tree.names=TRUE,digits=5)

  res
}

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



