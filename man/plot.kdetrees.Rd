\name{plot.kdetrees}
\alias{plot.kdetrees}
\title{Plot the unnormalized density estimates for each tree.}
\usage{
  \method{plot}{kdetrees} (x, ...)
}
\arguments{
  \item{x}{kdetrees object to be plotted}

  \item{...}{additional arguments passed to ggplot}
}
\value{
  a ggplot object
}
\description{
  Plot the unnormalized density estimates for each tree.
}
\examples{
fit <- kdetrees(apicomplexa)
plot(fit)
}
\author{
  Grady Weyenberg
}

