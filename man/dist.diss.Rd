\name{dist.diss}
\alias{dist.diss}
\title{pairwise tree distances}
\usage{
  dist.diss(x, ..., method = "euclidean", p = 2)
}
\arguments{
  \item{x}{either a row matrix of tree vectors, or a
  multiPhylo object}

  \item{...}{additional arguments passed to
  as.matrix.multiPhylo}

  \item{method}{option passed to dist}

  \item{p}{option passed to dist}
}
\value{
  a dist object with tree-to-tree distances
}
\description{
  pairwise tree distances
}
\examples{
dist.diss(apicomplexa[1:5])
}
\author{
  Grady Weyenberg
}
\seealso{
  dist
}

