\name{kdetrees.complete}
\alias{kdetrees.complete}
\title{Complete kdetrees analysis convenience function}
\usage{
  kdetrees.complete(file, ..., tree.file = "outliers.tre",
    csv.file = "results.csv", plot.file = "plot.png",
    hist.file = "hist.png")
}
\arguments{
  \item{file}{newick file with trees}

  \item{...}{additional parameters for kdetrees}

  \item{tree.file}{write outlier trees in newick format to
  this file}

  \item{csv.file}{write density results to this file}

  \item{plot.file}{print scatterplot of results to this
  file}

  \item{hist.file}{print histogram of density estimates to
  this file}
}
\value{
  result of kdetrees call
}
\description{
  Performs a complete kdetrees analysis, starting with
  reading trees from a newick file on disk, and writing
  result files to the working directory. Names and location
  of output files may be controlled by optional arguments.
}
\author{
  Grady Weyenberg
}

