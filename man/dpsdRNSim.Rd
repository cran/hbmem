\name{dpsdRNSim}
\alias{dpsdRNSim}
\title{Function dpsdRNSim}
\description{Simulate data from DPSD model with R a function of N}
\usage{
dpsdRNSim(NN = 2, NS = 1, I = 30, J = 200, K = 6, muN = c(-0.7, -0.5),
s2aN = 0.2, s2bN = 0.2, muS = 0, s2aS = 0.2, s2bS = 0.2,
muR = qnorm(0.25), phiA = -1, phiB = -1,
crit = matrix(rep(c(-1.6, -0.5, 0, 0.5, 1.6), each = I), ncol = (K - 1)))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{NN}{Number of new-item conditions.}
  \item{NS}{Number of studied-item conditions.}
  \item{I}{Number of participants.}
  \item{J}{Number of items.}
  \item{K}{Number of confidence ratings}
  \item{muN}{Mean of new-item distributuion}
  \item{s2aN}{Variance of participant effects on new-item distribution}
  \item{s2bN}{Variance of item effects on new-item distribution}
  \item{muS}{Mean of studied-item distribution}
  \item{s2aS}{Variance of participant effects on studied-item distribution}
  \item{s2bS}{Variance of item effects on studied-item distribution}
  \item{muR}{Mean of recollection (on probit space)}
  \item{phiA}{Linear slope of participant effect on recollection.}
  \item{phiB}{Linear slope of item effect on recollection.}
  \item{crit}{Matrix of criteria}
}
\references{See Pratte and Rouder (in review)}
