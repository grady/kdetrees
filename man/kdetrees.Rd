\name{kdetrees}
\alias{kdetrees}
\title{Identify discordant trees in a sample}
\usage{
  kdetrees(trees, n = ceiling(0.05 * length(trees)),
    bw = list(), ...)
}
\arguments{
  \item{trees}{multiPhylo object}

  \item{n}{number of outliers to report}

  \item{bw}{see Details}

  \item{...}{additional arguments for dist.diss}
}
\value{
  a kdetrees object; list(density,outliers,bandwidth)
}
\description{
  If bw is a single number, it will be used as a single
  constant bandwidth. It can also be a vector, in which
  case it will be used as variable bandwidths. Finally, if
  it is a list, the list will be passed as arguments to the
  bw.nn adaptive bandwith function.
}
\author{
  Grady Weyenberg
}

