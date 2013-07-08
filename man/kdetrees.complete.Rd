\name{kdetrees.complete}
\alias{kdetrees.complete}
\title{do complete analysis in one call}
\usage{
  kdetrees.complete(file, outgroup = NULL, ...,
    tree.file = "outliers.tre", csv.file = "results.csv",
    plot.file = "plot.png", hist.file = "hist.png")
}
\arguments{
  \item{file}{newick file with trees}

  \item{outgroup}{if a character, reroot all trees with
  this species as outgroup}

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
  results of kdetrees call
}
\description{
  do complete analysis in one call
}
\author{
  Grady Weyenberg
}

