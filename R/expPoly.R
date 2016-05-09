##' Pfaffian matrix for Exponential-polynomial distribution
##' 
##' This is the derivative matrix of the system. It should be arranged
##' as F_ij = dF(a)[j]/da[i]
##' @param th the vector theta of polynomial coeffiecients
##' @param G the value of the vector F(theta)
##' @return the derivative matrix evaluated at point theta
##' @author Grady Weyenberg
dG.expPoly <- function(th,G) {
    d <- length(th)
    b <- (1:d) * th
    
    x <- -(1+ b[-d]%*%G) / b[d]
    G <- c(G,x)
    M <- G[-1]
     for(i in 2:d){
        x <- -((i-1)*G[1] + b[-d] %*% G[-1])/b[d]
        G <- c(G[-1],x)
        M <- rbind(M,G[-1])
    }
    M
}



##' Evaluate 
##' 
##' Evaluate the higher order derivatives of an solved
##' Exponential-Polynomial system (Hayakawa and Takemura 2014).
##' @param th Parameter vector 
##' @param G solved F(th) vector
##' @param m order of derivative 
##' @return value
##' @author Grady Weyenberg
Gpartial <- function(th,G,m=0){
    ## m is derivitive to evaluate
    ## d is dimension of polynomial
    m <- trunc(m)
    if(m<0) stop("m must be nonnegative integer")
    d <- length(th)
    if(m < d - 1)
        res <- G[m+1]
    else {
        b <- (1:d)*th
        res <- -(1+b[-d]%*%G)/b[d]
        if(m >= d){
            G <- c(G,res)
            for(k in 2:(m-d+2)){
                res <- -((k-1)*G[1] + b[-d]%*%G[-1])/b[d]
                G <- c(G[-1],res)
            }
        }
    }
    res
}

##' Normalizing constant for Exponential-Polynomial distribution
##'
##' From Hayakawa and Takemura 2014
##' @param theta vector of polynomial coefficients
##' @param m derivative order of F to evaluate
##' @return Normalizing constant for the distribution
##' @author Grady Weyenberg
hgm.expPoly <- function(theta,m=0){
    if(all(theta>0))
        theta <- -theta
    if(any(theta>0))
        stop("all elements of theta must be negative")
    d <- length(theta)
    th0 <- rep(0.0,d)
    th0[d] <- -1.0 ##theta[d]
    m <- (1:(d-1)) / d
    G0 <- gamma(m)/d #/ (-th0[d])^m #Eqn (9) Hayakawa & Takemura 2014
    G1 <- hgm.Rhgm(th0,G0,theta,dG.expPoly)
    ## Gpartial(theta,G1,m)
    G1
}

const.expPoly <- function(theta,m=0){
    exp.poly <- function(x,theta){ ##Eqn (5) from paper in review.
        y <- 1.0
        ## res <- x
        ## res[TRUE] <- 0.0
        res <- numeric(length(x)) #zeros
        for(t in theta){
            y <- y * x
            res <- res + t * y
        }
        x^m * exp(res)
    }
    G1 <- integrate(exp.poly,0,Inf,theta=theta)$value
    ##Gpartial(theta,G1,m) ##This will only work with length(theta)==2
    G1
}

##' Create a theta vector from a phylo object
##'
##' Converts a tree into a theta vector for the exponential polynomial
##' distribution needed for normalizing constant calculation.
##' @param t tree (a phylo object)
##' @param h bandwidth
##' @return theta vector for expPoly distribution
##' @author Grady Weyenberg
make.theta <- function(t,h){
    p <- length(t$edge.length)
    d0t.sq <- crossprod(t$edge.length)
    -c(2.0*sqrt(d0t.sq), 1.0) / h^2
}

