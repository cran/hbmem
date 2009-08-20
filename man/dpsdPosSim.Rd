\name{dpsdPosSim}
\alias{dpsdPosSim}
\title{Function dpsdPosSim}
\description{Function to simulate data from DPSD model with positive d'}
\usage{
dpsdPosSim(I = 30, J = 200, K = 6, muN = -0.7, s2aN = 0.2, s2bN = 0.2, muD = 0, s2aD = 0.2, s2bD = 0.2, muR = qnorm(0.25), s2aR = 0.2, s2bR = 0.2, crit = matrix(rep(c(-1.6, -0.5, 0, 0.5, 1.6), each = I), ncol = (K - 1)))
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
  \item{s2aR}{ ~~Describe \code{s2aR} here~~ }
  \item{s2bR}{ ~~Describe \code{s2bR} here~~ }
  \item{crit}{ ~~Describe \code{crit} here~~ }
}
\keyword{internal}