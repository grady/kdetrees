##' Extract vector of internal branch lengths
##'
##' Take a phylo object and extract a vector containing the internal
##' edge lengths, sorted by magnitude.
##' @param x a phylo object
##' @return a vector with the internal edge lengths
##' @author Grady Weyenberg
internal.blens <- function(x){
    x <- unroot(x)
    n <- Ntip(x)
    has.leaf.node <- function(z) {any(z %in% seq_len(n))}
    i <- apply(x$edge,1,has.leaf.node)
    sort(x$edge.length[!i])
}

##' Count semi-labeled rooted trees on n taxa
##'
##' Counts the number of topologies for rooted trees on n taxa. This
##' is the number of orthants in BHV tree space.
##' 
##' @param n number of taxa, must be at least 3
##' @return integer number of topologies for a rooted phylogeny
##' @author Grady Weyenberg
Northant <- function(n){
  n <- trunc(n)
  if (any(n<3)) stop("n must be greater than 2")
  x <- 2 * n - 3
  res <- rep(1,length(n))
  while (any(x > 1)) {  # double factorial
    res <- res * x
    x <- x - 2
    x[x<1] <- 1
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
bhv.consts <- function(trees,bw, hgm=TRUE){
  x <- sapply(trees,internal.blens)
  if(!is.matrix(x)) stop("different number of tips in trees")
  n <- Northant(Ntip(trees[[1]])) # number of orthants
  d <- NROW(x) # dim of orthants
  N <- NCOL(x) # number of trees
  res <- numeric(N) # result obj
  
  if ( length(bw) == 1L ) bw <- rep(bw,N)
  
  for (i in seq_len(N)){
    bm <- diag(d) * bw[i] * bw[i]
    c0 <- 1/sqrt((2 * pi)^d * det(bm)) # gaussian const at center
    if(min(x[,i]) > 6.0 * bw[i]){
        res[i] <- c0
    } else {
        p0 <- bhv.orthant.lb(trees[[i]], bw[i])
        if(hgm){
            p1 <- hgm.ncorthant(bm, x[,i]) / c0
            res[i] <- (p1 + p0*(n-1))
        } else {
            res[i] <- p0 * n
        }
    }
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


##' Lower orthant bound for tree t with bw h.
##' @param t tree
##' @param h bandwidth
##' @return lower bound volmue
##' @author Grady Weyenberg
bhv.orthant.lb <- function(t,h){
    e <- internal.blens(t)
    p <- length(e)
    d0t.sq <- crossprod(e)
    theta <- -c(2.0*sqrt(d0t.sq), 1.0) / h^2
    A <- const.expPoly(theta,p-1)
    pi^(p/2) * exp(-d0t.sq/h^2) * A / 2^(p-1) / gamma(p/2)
}
