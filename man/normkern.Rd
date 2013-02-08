\name{normkern}
\alias{normkern}
\title{Generalized Gaussian kernel}
\usage{
  normkern(x, bw = 1, delta = 2L)
}
\arguments{
  \item{x}{places to evaluate kernel}

  \item{bw}{bandwidth values}

  \item{delta}{shape parameter for kernel}
}
\value{
  an object of the same type as x with the kernel
  evaluations
}
\description{
  The un-normalized Gaussian kernel function:
  exp(-(abs(x/bw))^delta)/bw
}
\author{
  Grady Weyenberg
}

