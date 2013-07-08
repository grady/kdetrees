\name{as.matrix.multiPhylo}
\alias{as.matrix.multiPhylo}
\title{dissimilarity map tree vectorization}
\usage{
  \method{as.matrix}{multiPhylo} (x, use.blen = FALSE,
    unroot = TRUE, ...)
}
\arguments{
  \item{x}{an ape::multiPhylo object.}

  \item{use.blen}{use branch length information in the
  calculation?}

  \item{unroot}{unroot trees first?}

  \item{...}{additional options for \code{dist}}
}
\value{
  a row matrix of tree vectors
}
\description{
  Dissimilarity maps convert trees to vectors using
  tip-to-tip path lengths. Branch length information may be
  optionally discarded (the default), resulting in vectors
  based solely on tree topology.
}
\author{
  Grady Weyenberg
}

