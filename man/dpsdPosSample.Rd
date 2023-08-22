\name{dpsdPosSample}
\alias{dpsdPosSample}
\title{Function dpsdPosSample}
\description{Function to sample posterior of DPSD model with positive d'}
\usage{
dpsdPosSample(dat, M = 5000, keep = (M/10):M, getDIC = TRUE, jump = 0.01)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{Data frame that must include variables
    Scond,cond,sub,item,lag,resp.  Scond indexes studied/new, whereas
    cond indexes conditions nested within the studied or new
    conditions.  Indexes for Scond,cond, sub, item, and respone must
    start at zero and have no gaps (i.e., no skipped subject numbers). Lags
    must be zero-centered.}
  \item{M}{Number of MCMC iterations.}
  \item{keep}{Which MCMC iterations should be included in estimates and
    returned.  Use keep to both get ride of burn-in, and thin chains if
    necessary.}
  \item{getDIC}{Logical. Should the function compute DIC value?  This
    takes a while if M is large.}
  \item{jump}{The criteria and decorrelating steps utilize
    Matropolis-Hastings sampling routines, which require tuning.  All
    MCMC functions should self-tune during the burnin period (iterations
    before keep), and they will alert you to the success of tuning.  If
    acceptance rates are too low, "jump" should be decreased, if they
    are too hight, "jump" should be increased.  Alternatively, or in
    addition to adjusting "jump", simply increase the burnin period
    which will allow the function more time to self-tune.}
}
\keyword{internal}