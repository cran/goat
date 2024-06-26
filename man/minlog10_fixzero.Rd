% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/util_stats.R
\name{minlog10_fixzero}
\alias{minlog10_fixzero}
\title{-log10 transform a vector of p-values, replacing zeros with some limit/threshold}
\usage{
minlog10_fixzero(x, limit = 2.22e-16)
}
\arguments{
\item{x}{p-value vector to transform to -log10}

\item{limit}{value to replace zero's in \code{x} with. Set NA to replace zero's in \code{x} with the smallest finite value in \code{x} (if there is none, defaults to 2.22e-16)}
}
\value{
input parameter \code{x} transformed to -log10
}
\description{
-log10 transform a vector of p-values, replacing zeros with some limit/threshold
}
\examples{
  pval = c(0, 10^-6, 0.001, 0.01, 1, NA, -Inf, Inf, NaN)
  cbind(
    input = pval,
    # default; replace zeros with typical R machine precision for doubles
    minlog10_default = minlog10_fixzero(pval),
    # alternatively, replace zero with lowest non-zero pvalue in input
    minlog10_limit_from_data = minlog10_fixzero(pval, limit = NA)
  )
}
