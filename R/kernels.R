### Copyright (C) 2014 -- Grady Weyenberg ###
### Kernel functions are defined here. These functions have a template
### function(x,bw,...)
### they should be symmetric about zero and should integrate to a constant
### value for all positive bw.

##' Generalized Gaussian kernel
##'
##' The un-normalized Gaussian kernel function: exp(-(abs(x/bw))^delta)/bw
##'
##' The bandwidth parameter may be used in any way that makes sense in
##' the above R expression. In particular, it may be a single value,
##' for a constant bandwidth, or a vector, with each element
##' corresponding the bandwidth of the kernel to be placed at each
##' respective observation.
##' 
##' @param x places to evaluate kernel
##' @param bw bandwidth values
##' @param bhv.c vector of kernel constants 
##' @param delta shape parameter for kernel
##' @return an object of the same type as x with the kernel evaluations
##' @author Grady Weyenberg
normkern <- function(x, bw=1.0, bhv.c=NULL, delta=2L){
    if(is.null(bhv.c)){
      exp(-abs(x/bw)^delta) / bw
    } else {
      exp(-abs(x/bw)^delta) / bhv.c
    }
}
    

## place a kernel of given bw on each observation (rows?cols?)
## each kernel neds to be normalized, in geodesic or regular Rn
## in geodesic space, in geodesic space it will be

foo <- function(tree,bw){
    bl <- internal.blens(tree)
    p1 <- hgm.ncorthant(diag(length(bl)) * bw, bl)
    p0 <- hgm.ncorthant(diag(length(bl)) * bw, -bl)
}
