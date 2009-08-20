\name{gammaLikeSim}
\alias{gammaLikeSim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Function gammaLikeSim}
\description{Simulates data from gamma model with criteria on likelihood
space.}
\usage{
gammaLikeSim(NN = 1, NS = 1, I = 30, J = 200, K = 6, muS = log(0.5), s2aS = 0.2, s2bS = 0.2, lagEffect = 0, shape = 2, crit = matrix(rep(c(0.75, 0.8, 1, 1.35, 1.6), each = I), ncol = (K - 1)))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{NN}{ ~~Describe \code{NN} here~~ }
  \item{NS}{ ~~Describe \code{NS} here~~ }
  \item{I}{ ~~Describe \code{I} here~~ }
  \item{J}{ ~~Describe \code{J} here~~ }
  \item{K}{ ~~Describe \code{K} here~~ }
  \item{muS}{ ~~Describe \code{muS} here~~ }
  \item{s2aS}{ ~~Describe \code{s2aS} here~~ }
  \item{s2bS}{ ~~Describe \code{s2bS} here~~ }
  \item{lagEffect}{ ~~Describe \code{lagEffect} here~~ }
  \item{shape}{ ~~Describe \code{shape} here~~ }
  \item{crit}{ ~~Describe \code{crit} here~~ }
}
\keyword{internal}