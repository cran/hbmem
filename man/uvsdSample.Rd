\name{uvsdSample}
\alias{uvsdSample}
\title{Function uvsdSample}
\description{Runs MCMC for the hierarchical UVSD model}
\usage{uvsdSample(dat, M = 10000, keep = (M/10):M, getDIC = TRUE,
freeCrit=TRUE, equalVar=FALSE, freeSig2=FALSE, jump=.01)}

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
  \item{freeCrit}{Logical.  If TRUE (default) individual criteria vary
    across people.  If false, all participants have the same criteria
    (but note that overall response biases are still modeled in the
    means)}
  \item{equalVar}{Logical. If FALSE (default), unequal-variance model
  is fit.  If TRUE, equal-variance model is fit.}
  \item{freeSig2}{Logical.  If FALSE (default), one sigma is fit
  for all participants and items (as in Pratte, et al., 2009).  If TRUE,
  then an additive model is placed on the log of sigma2.}
  \item{jump}{The criteria and decorrelating steps utilize
    Matropolis-Hastings sampling routines, which require tuning.  All
    MCMC functions should self tune during the burnin perior (iterations
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
  \item{estS2}{Posterior means of block for studied-item variances.}
  \item{estCrit}{Posterior means of criteria}
  \item{blockN}{Each iteration for each parameter in the new-item mean
  block.  Rows index iteration, columns index parameter.}
  \item{blockS}{Same as blockN, but for the studied-item means}
  \item{blockS2}{Same as blockN, but for variances of studied-item
  distribution.  If equalVar=TRUE, then these values are all zero.  If
  UVSD is fit but freeSig2=FALSE, then only the first element is
  non-zero (mu).}
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
    keep, or jump need to be adjusted.}
   \item{b0S2}{If additive model is placed on Sigma2 (i.e.,
   freeSigma2=TRUE),  then all parameters on S2 must be tuned.  b0S2 are
   the acceptance probabilities for these parameters.}
}
\references{See Pratte, Rouder, & Morey (2009)}
\author{Michael S. Pratte}
\seealso{hbmem}
\examples{
#In this example we generate data from EVSD, then fit it with both
#hierarchical UVSD and UVSD assuming no participant or item effects.
#Note that there is variance in the studied-item distribution, but not
#the new-item distribution.

library(hbmem)
sim=uvsdSim(I=70,J=400,s2aN = 0, s2bN = 0,s2aS=.25,s2bS=.25)
dat=as.data.frame(cbind(sim@subj,sim@item,sim@cond,sim@lag,sim@resp))
colnames(dat)=c("sub","item","cond","lag","resp")

M=100
UVSD=uvsdSample(dat,M=M)
#Now make data frame for data-averaged analysis
adat=dat
adat$sub=0
adat$item=0
aUVSD=uvsdSample(adat,M=M,freeCrit=FALSE)

par(mfrow=c(1,2))
plot(sqrt(exp(UVSD@blockS2[,UVSD@mu])),t='l',main="Hierarchical")
plot(sqrt(exp(aUVSD@blockS2[,aUVSD@mu])),t='l',main="Data-Averaged")
}
