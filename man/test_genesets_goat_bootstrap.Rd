% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/test_genesets_goat.R
\name{test_genesets_goat_bootstrap}
\alias{test_genesets_goat_bootstrap}
\title{Naive GOAT variant where we estimate null parameters for each geneset size independently}
\usage{
test_genesets_goat_bootstrap(
  genesets,
  genelist,
  score_type,
  niter = 5e+05,
  verbose = FALSE
)
}
\arguments{
\item{genesets}{see \code{test_genesets_goat_precomputed}}

\item{genelist}{see \code{test_genesets_goat_precomputed}}

\item{score_type}{see \code{test_genesets_goat_precomputed}}

\item{niter}{integer number of bootstrap iterations; at least 10000, at most 5000000}

\item{verbose}{boolean, create debug plots}
}
\value{
see \code{test_genesets_goat_precomputed}
}
\description{
In typical use-cases, one applies \code{test_genesets()} instead with parameter \code{method="goat"} ,
which in turn will use \code{test_genesets_goat_precomputed} (and not this function).
}
