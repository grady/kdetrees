### Kernel functions are defined here. These functions have a template
### function(x,bw,...)
### they should be symmetric about zero and should integrate to a constant
### value for all positive bw.

##' Generalized Gaussian kernel
##'
##' The un-normalized Gaussian kernel function: exp(-(abs(x/bw))^delta)/bw
##' 
##' @param x places to evaluate kernel
##' @param bw bandwidth values
##' @param delta shape parameter for kernel
##' @return an object of the same type as x with the kernel evaluations
##' @author Grady Weyenberg
normkern <- function(x, bw=1.0, delta=2L)
  exp(-abs(x/bw)^delta) / bw

