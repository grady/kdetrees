% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/expPoly.R
\name{make.theta}
\alias{make.theta}
\title{Create a theta vector from a phylo object}
\usage{
make.theta(t, h)
}
\arguments{
\item{t}{tree (a phylo object)}

\item{h}{bandwidth}
}
\value{
theta vector for expPoly distribution
}
\description{
Converts a tree into a theta vector for the exponential polynomial
distribution needed for normalizing constant calculation.
}
\author{
Grady Weyenberg
}

