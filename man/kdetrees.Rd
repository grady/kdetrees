\name{kdetrees}
\alias{kdetrees}
\title{Identify discordant trees in a sample}
\usage{
  kdetrees(trees, k = 1.5,
    distance = c("geodesic", "dissimilarity"), bw = list(),
    greedy = FALSE, ...)
}
\arguments{
  \item{trees}{multiPhylo object}

  \item{k}{IQR multiplier for outlier detection}

  \item{distance}{Select "geodesic" or "dissimilarity"
  distance calculation method}

  \item{bw}{see Details}

  \item{greedy}{greedy outlier detection?}

  \item{...}{additional arguments for distance calculation
  function, see details}
}
\value{
  a kdetrees object; list(density,outliers)
}
\description{
  If bw is a single number, it will be used as a single
  constant bandwidth. It can also be a vector, in which
  case it will be used as variable bandwidths. Finally, if
  it is a list, the list will be passed as arguments to the
  bw.nn adaptive bandwith function.
}
\details{
  ... Is passed to either \code{distory::dist.multiPhylo}
  or \code{dist.diss}, as appropriate. See the help for
  these functions for more details.
}
\examples{
kdeobj <- kdetrees(apicomplexa)
print(kdeobj)
kdeobj$outliers

kdetrees(apicomplexa, k=2.0, distance="dissimilarity",use.blen=TRUE)
}
\author{
  Grady Weyenberg
}

