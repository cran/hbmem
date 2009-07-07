\name{dpsdSample}
\alias{dpsdSample}
\title{Function dpsdSample}
\description{Runs MCMC for the hierarchical DPSD model}
\usage{dpsdSample(dat, M = 5000, keep = (M/10):M, getDIC = TRUE, jump=.01)}

\arguments{
  \item{dat}{Data frame that must include variables
    cond,sub,item,lag,resp.  Indexes for cond, sub, item, and respone
    must start at zero and have no gapes (i.e., no skipped subject
    numbers). Lags must be zero-centered.}
  \item{M}{Number of MCMC iterations.}
  \item{keep}{Which MCMC iterations should be included in estimates and
    returned.  Use keep to both get ride of burn-in, and thin chains if necessary}
  \item{getDIC}{Logical. should the function compute DIC value?  This
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
\value{
  The function returns an internally defined "uvsd" structure that
  includes the following components
  \item{mu}{Indexes which element of blocks contain mu}
  \item{alpha}{Indexes which element of blocks contain participant
    effects, alpha}
  \item{beta}{Indexes which element of blocks contain item effects, beta}
  \item{s2alpha}{Indexes which element of blocks contain variance of
    participant effects (alpha).}
  \item{s2beta}{Indexes which element of blocks contain variance of
    item effects (beta).}
  \item{theta}{Indexes which element of blocks contain theta, the slope of
    the lag effect}
  \item{estN}{Posterior means of block parameters for new-item means}
  \item{estS}{Posterior means of block parameters for studied-item means}
  \item{estR}{Posterior means of block for Recollection means.}
  \item{estCrit}{Posterior means of criteria}
  \item{blockN}{Each iteration for each parameter in the new-item mean
  block.  Rows index iteration, columns index parameter.}
  \item{blockS}{Same as blockN, but for the studied-item means}
  \item{blockR}{Same as blockN, but for the recollection-parameter
    means.}
  \item{s.crit}{Samples of each criteria.}
  \item{pD}{Number of effective parameters used in DIC.  Note that this
    should be smaller than the actual number of parameters, as
    constraint from the hierarchical structure decreases the number of
    effective parameters.}
  \item{DIC}{DIC value.  Smaller values indicate better fits.  Note that
    DIC is notably biased toward complexity.}
  \item{M}{Number of MCMC iterations run}
  \item{keep}{MCMC iterations that were used for estimation and
    returned}
  \item{b0}{Metropolis-Hastings acceptance rates for decorrelating
    steps.  These should be between .2 and .6.  If they are not, the M,
    keep, or jump arguments need to be adjusted.}
  \item{b0Crit}{acceptance rates for criteria.}
}
\references{See Pratte, Rouder, & Morey (2009)}
\author{Michael S. Pratte}
\seealso{hbmem}
\examples{
#In this example we generate data from EVSD, then fit it with both
#hierarchical DPSD and DPSD assuming no participant or item effects.
library(hbmem)
sim=dpsdSim(I=50,J=200)
dat=as.data.frame(cbind(sim@subj,sim@item,sim@cond,sim@lag,sim@resp))
colnames(dat)=c("sub","item","cond","lag","resp")

M=500 #Too low for real analysis!
DPSD=dpsdSample(dat,M=M)

#Look at recollection estimates
par(mfrow=c(1,3),pch=19,pty='s')
plot(pnorm(DPSD@blockR[,DPSD@mu]),t='l',
ylab="P(recollection)")
abline(h=pnorm(sim@muR),col="blue")
plot(DPSD@estR[DPSD@alpha]~sim@alphaR)
abline(0,1,col="blue")
plot(DPSD@estR[DPSD@beta]~sim@betaR)
abline(0,1,col="blue")
}
