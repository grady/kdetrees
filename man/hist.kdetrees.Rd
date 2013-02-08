\name{hist.kdetrees}
\alias{hist.kdetrees}
\title{Create a histogram of tree density estimates}
\usage{
  \method{hist}{kdetrees} (x, ...,
    ggplot = ("ggplot2" \%in\% installed.packages()))
}
\arguments{
  \item{x}{kdetrees object to plot}

  \item{...}{see details}

  \item{ggplot}{Logical: ggplot2 or base graphics}
}
\value{
  either a ggplot or histogram object
}
\description{
  If ggplot2 graphics are used, ... is a sequence of
  additional ggplot2 directives to be added to the end of
  the ggplot construction. If base graphics are used, ...
  are additional arguments passed to plot command.
}
\examples{
fit <- kdetrees(apicomplexa,12,use.blen=TRUE)
plot(fit)
}
\author{
  Grady Weyenberg
}

