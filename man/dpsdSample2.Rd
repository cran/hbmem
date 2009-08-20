\name{dpsdSample2}
\alias{dpsdSample2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Fit DPSD model with means in participant effect}
\description{
This is just a way to sample the DPSD model with no mus (grand means) to
make it easier to fit the models to data with no item effects and
non-hierarchical priors.  In short, this should never be used except for
comparison purposes.
}
\usage{
dpsdSample2(dat, M = 5000, keep = (M/10):M, getDIC = TRUE, Hier = TRUE, jump = 0.01)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{ ~~Describe \code{dat} here~~ }
  \item{M}{ ~~Describe \code{M} here~~ }
  \item{keep}{ ~~Describe \code{keep} here~~ }
  \item{getDIC}{ ~~Describe \code{getDIC} here~~ }
  \item{Hier}{ ~~Describe \code{Hier} here~~ }
  \item{jump}{ ~~Describe \code{jump} here~~ }
}
\keyword{internal}