\name{dist.diss}
\alias{dist.diss}
\title{pairwise tree distances}
\usage{
  dist.diss(x, ..., method = "euclidean", p = 2)
}
\arguments{
  \item{x}{either a row matrix of tree vectors, or a
  multiPhylo object}

  \item{...}{additional arguments passed to dissmap}

  \item{method}{option passed to dist}

  \item{p}{option passed to dist}
}
\value{
  a matrix of pairwise tree distances
}
\description{
  pairwise tree distances
}
\author{
  Grady Weyenberg
}
\seealso{
  dist
}

