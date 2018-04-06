\name{dpsdFullSim}
\alias{dpsdFullSim}
\title{Function dpsdFullSim}
\description{Simulate data from DPSD model with R a function of d' plus
  noise.  Model not yet published.}
\usage{
dpsdFullSim(I = 30, J = 200, K = 6, muN = -0.7, s2aN = 0.2,
s2bN = 0.2, muD = 0, s2aD = 0.2, s2bD = 0.2, muR = qnorm(0.25),
phiA = 1, etaA = 0, phiB = 1, etaB = 0,crit =
matrix(rep(c(-1.6, -0.5, 0, 0.5, 1.6), each = I), ncol = (K - 1)))
}
\arguments{
  \item{I}{ ~~Describe \code{I} here~~ }
  \item{J}{ ~~Describe \code{J} here~~ }
  \item{K}{ ~~Describe \code{K} here~~ }
  \item{muN}{ ~~Describe \code{muN} here~~ }
  \item{s2aN}{ ~~Describe \code{s2aN} here~~ }
  \item{s2bN}{ ~~Describe \code{s2bN} here~~ }
  \item{muD}{ ~~Describe \code{muD} here~~ }
  \item{s2aD}{ ~~Describe \code{s2aD} here~~ }
  \item{s2bD}{ ~~Describe \code{s2bD} here~~ }
  \item{muR}{ ~~Describe \code{muR} here~~ }
  \item{phiA}{ ~~Describe \code{phiA} here~~ }
  \item{etaA}{ ~~Describe \code{etaA} here~~ }
  \item{phiB}{ ~~Describe \code{phiB} here~~ }
  \item{etaB}{ ~~Describe \code{etaB} here~~ }
  \item{crit}{ ~~Describe \code{crit} here~~ }
}
\keyword{internal}