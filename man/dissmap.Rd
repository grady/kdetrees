\name{dissmap}
\alias{dissmap}
\title{dissimilarity map tree vectorization}
\usage{
  dissmap(x, use.blen = FALSE, unroot = TRUE, ...)
}
\arguments{
  \item{x}{an ape::multiPhylo object. All trees must have
  identical tip names.}

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

