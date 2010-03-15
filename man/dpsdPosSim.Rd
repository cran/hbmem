\name{dpsdPosSim}
\alias{dpsdPosSim}
\title{Function dpsdPosSim}
\description{Function to simulate data from DPSD model with positive d'}
\usage{
dpsdPosSim(NN = 1, NS = 2, I = 30, J = 200, K = 6, muN = -0.7,
                 s2aN = 0.2, s2bN = 0.2, muD = c(0, 0.5), s2aD = 0.2,
                 s2bD = 0.2, muR = qnorm(c(0.2, 0.4)), s2aR = 0.2, s2bR
                 = 0.2, crit = matrix(rep(c(-1.6, -0.5, 0, 0.5, 1.6),
                 each = I), ncol = (K - 1)))
}
\arguments{
  \item{NN}{Number of new-item conditions.}
  \item{NS}{Number of studied-item conditions.}
  \item{I}{Number of participants.}
  \item{J}{Number of items.}
  \item{K}{Number of response options.}
  \item{muN}{Mean of new-item distribution.  If there are more than one new-item conditions this is a vector of means with length equal to NN.}
  \item{s2aN}{Variance of participant effects on mean of new-item distribution.}  
  \item{s2bN}{Variance of item effects on mean of new-item distribution.}
  \item{muD}{Mean of dprime distribution.  If there are more than
    new-item conditions this is a vector of means with length equal to NNone studied-item conditions this is a vector of means with length equal to NS.}
  \item{s2aD}{Variance of participant effects on mean of dprime distribution.}  
  \item{s2bD}{Variance of item effects on mean of dprime distribution.}
  \item{muR}{Mean recollection, on probit space.}
  \item{s2aR}{Variance of participant effects recollection.}
  \item{s2bR}{Variance of item effects on recollection.}
  \item{crit}{Matrix of criteria (not including -Inf or Inf).  Columns
    correspond to criteria, rows correspond to participants.}
}
\keyword{internal}