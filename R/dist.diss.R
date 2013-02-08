##' dissimilarity map tree vectorization
##'
##' Dissimilarity maps convert trees to vectors using tip-to-tip path
##' lengths. Branch length information may be optionally discarded
##' (the default), resulting in vectors based solely on tree
##' topology.
##' 
##' @param x an ape::multiPhylo object. All trees must have identical
##' tip names.
##' @param use.blen use branch length information in the calculation?
##' @param unroot unroot trees first?
##' @param ... additional options for \code{dist}
##' @return a row matrix of tree vectors
##' @author Grady Weyenberg
##' @export
dissmap <- function(x,use.blen=FALSE,unroot=TRUE,...){
  ##Check if all tres have same tips
  tip.labels <- sort(x[[1]]$tip.label)
  stopifnot( tip.labels == sapply(x,function(y) sort(y$tip.label)) )
  ##Unroot the trees (is there any reason not to do this?)
  if(unroot){ x <- lapply(x,unroot) }
  ##Set branch lengths to 1
  if(!use.blen){ x <- lapply(x, compute.brlen, method = 1) }
  ##find tip-tip distances for a tree
  tip2tip <- function(y){
    o <- cophenetic(y,...)[tip.labels,tip.labels]
    o[upper.tri(o)]
  }
  cnames <- outer(tip.labels,tip.labels,paste,sep="-")
  ##convert multiPhylo to a row matrix of tip distances
  out <- t(sapply(x,tip2tip))
  colnames(out) <- cnames[upper.tri(cnames)]
  if(is.null(rownames(out))) rownames(out) <- paste("tree",1:nrow(out),sep="")
  out
}

##' pairwise tree distances
##'
##' @param x either a row matrix of tree vectors, or a multiPhylo object
##' @param ... additional arguments passed to dissmap
##' @param method option passed to dist
##' @param p option passed to dist
##' @seealso dist
##' @return a matrix of pairwise tree distances
##' @author Grady Weyenberg
##' @export
dist.diss <- function(x,...,method="euclidean",p=2){
  ##pairwise stree distances
  if(inherits(x,"multiPhylo")) d <- dissmap(x,...)
  as.matrix(dist(d,method=method,p=p))
} 
