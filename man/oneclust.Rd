% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/oneclust.R
\name{oneclust}
\alias{oneclust}
\title{Maximum homogeneity clustering for one-dimensional data}
\usage{
oneclust(x, k, w = NULL, sort = TRUE)
}
\arguments{
\item{x}{Numeric vector, samples to be clustered.}

\item{k}{Integer, number of clusters.}

\item{w}{Numeric vector, sample weights (optional).
Note that the weights here should be sampling weights
(for example, a certain proportion of the population),
not frequency weights (for example, number of occurrences).}

\item{sort}{Should we sort \code{x} (and \code{w}) before clustering?
Default is \code{TRUE}. Otherwise the order of the data is respected.}
}
\value{
A list containing:
\itemize{
\item \code{cluster} - cluster id of each sample.
\item \code{cut} - index of the optimal cut points.
}
}
\description{
Maximum homogeneity clustering for one-dimensional data
}
\examples{
set.seed(42)
x <- sample(c(
  rnorm(50, sd = 0.2),
  rnorm(50, mean = 1, sd = 0.3),
  rnorm(100, mean = -1, sd = 0.25)
))
oneclust(x, 3)
}
\references{
Fisher, Walter D. 1958. On Grouping for Maximum Homogeneity.
\emph{Journal of the American Statistical Association} 53 (284): 789--98.
}
