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

##' Count semi-labeled rooted trees on n taxa
##'
##' Counts the number of topologies for rooted trees on n taxa. This
##' is the number of orthants in BHV tree space.
##' (NOTE) Candidate for optimization?
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

##' Find BHV constant lower bound
##'
##' Obtain a lower bound for the kernel integration constants
##' for centered on the trees with bandwidths bw.
##' @param trees list of trees at which kernels will be centered
##' @param bw bandwidth for each kernel
##' @return a vector of kernel integration constant lower bounds
##' @author Grady Weyenberg
bhv.consts <- function(trees,bw){
  ## x <- sapply(trees,internal.blens) ##dist.multiPhylo uses terminal edges
  x <- sapply(trees,getElement,"edge.length")
  if(!is.matrix(x)) stop ("different number of tips in trees")
  n <- Northant(Ntip(trees[[1]])) #number of orthants
  d <- NROW(x) #dim of orthants
  N <- NCOL(x) #number of trees
  res <- numeric(N) #result obj
  
  if ( length(bw) == 1L ) bw <- rep(bw,N)
  
  for (i in 1:N){
    bm <- diag(d) * bw[i] * bw[i]
    c0 <- 1/sqrt((2 * pi)^d * det(bm)) #gaussian const at center
    ##browser()
    ##p1 <- hgm.ncorthant(bm, x[,i]) / c0
    p0 <- bhv.orthant.lb(trees[[i]], bw[i])
    ##res[i] <- (p1 + p0*(n-1))
    res[i] <- p0*n
  }
  res
}

## bhv.consts.ub <- function(trees,bw){
##     x <- sapply(trees,internal.blens)
##     if(!is.matrix(x)) stop ("different Ntips in trees")
##     n <- Northant(Ntip(trees[[1]])) #number of orthants
##     d <- nrow(x) #dim of orthants
##     N <- ncol(x) #number of trees
##     res <- numeric(N) #result obj

##     if ( length(bw) == 1L ) bw <- rep(bw,N)

##     for (i in 1:N){
##         bm <- diag(d) * bw[i]
##         c0 <- 1/sqrt((2 * pi)^d * det(bm)) #gaussian const at center
##         p1 <- hgm.ncorthant(bm, x[,i])
##         p0 <- hgm.ncorthant(bm, -x[,1])
##         y <- x[,i]
##         y[1] <- -y[1]
##         p2 <- hgm.ncorthant(bm,y)
        
##         res[i] <- (p1 + p2 * 2 + p0 * (n-3))/c0
##     }
##     res
## }