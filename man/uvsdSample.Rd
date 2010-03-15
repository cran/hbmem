\name{uvsdSample}
\alias{uvsdSample}
\title{Function uvsdSample}
\description{Runs MCMC estimation for the hierarchical UVSD model.}
\usage{uvsdSample(dat, M = 10000, keep = (M/10):M, getDIC = TRUE,
freeCrit=TRUE, equalVar=FALSE, freeSig2=FALSE, Hier=TRUE,jump=.0001)}

\arguments{
  \item{dat}{Data frame that must include variables
    Scond,cond,sub,item,lag,resp.  Scond indexes studied/new, whereas
    cond indexes conditions nested within the studied or new
    conditions.  Indexes for Scond,cond, sub, item, and response must
    start at zero and have no gaps (i.e., no skipped subject numbers). Lags
    must be zero-centered.}
  \item{M}{Number of MCMC iterations.}
  \item{keep}{Which MCMC iterations should be included in estimates and
    returned.  Use keep to both get ride of burn-in, and thin chains if necessary}
  \item{getDIC}{Logical. should the function compute DIC value?  This
    takes a while if M is large.}
  \item{freeCrit}{Logical.  If TRUE (default) individual criteria vary
    across people.  If false, all participants have the same criteria.
    This should be set to false if there is only one participant, e.g.,
    if averaging data over subjects.}
  \item{equalVar}{Logical. If FALSE (default), unequal-variance model
  is fit.  If TRUE, equal-variance model is fit.}
  \item{freeSig2}{Logical.  If FALSE (default), one sigma is fit
  for all participants and items (as in Pratte, et al., 2009).  If TRUE,
  then an additive model is placed on the log of sigma2 (as in Pratte
  and Rouder (2010).}
\item{Hier}{Logical. If TRUE then the variances of effects
    (e.g., item effects) are estimated from the data, i.e., effects are
    treated as random.  If FALSE then these variances are fixed to
    2.0 (.5 for recollection effects), thus treating these effects as
    fixed.  This option is there to allow for compairson with more
    traditional approaches, and to see the effects of imposing
    hierarcical structure.  It should always be set to TRUE in real
    analysis, and is not even guaranteed to work if set to false.}
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
  The function returns an internally defined "uvsd" S4 class that
  includes the following components
  \item{mu}{Indexes which element of blocks contain grand means, mu}
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
#In this example we generate data from UVSD with a different muN,muS,and
#Sigma2 for every person and item. These data are then fit with 
#hierarchical UVSD allowing participant or item effects on log(sigma2).

library(hbmem)
sim=uvsdSim(NN=1,muN=-.5,NS=2,muS=c(.5,1),I=30,J=300,s2aN = .2, s2bN = .2,
muS2=log(c(1.3,1.5)),s2aS=.2,s2bS=.2,s2aS2=.2,s2bS2=.2)
dat=as.data.frame(cbind(sim@subj,sim@item,sim@cond,sim@Scond,sim@lag,sim@resp))
colnames(dat)=c("sub","item","cond","Scond","lag","resp")

M=10 #Way too low for real analysis
keep=2:M
uvsd=uvsdSample(dat,M=M,keep=keep,equalVar=FALSE,freeSig2=TRUE,jump=.0001,Hier=1)

par(mfrow=c(3,2),pch=19,pty='s')
#Look at chains of MuN and MuS
matplot(uvsd@blockN[,uvsd@muN],t='l',xlab="Iteration",ylab="Mu-N")
abline(h=sim@muN,col="blue")
matplot(uvsd@blockS[,uvsd@muS],t='l',xlab="Iteration",ylab="Mu-S")
abline(h=sim@muS,col="blue")

#Estimates of strength effects as function of true values
plot(uvsd@estN[uvsd@alphaN]~sim@alphaN,xlab="True
Alpha-N",ylab="Est. Alpha-N");abline(0,1,col="blue")
plot(uvsd@estS[uvsd@alphaS]~sim@alphaS,xlab="True
Alpha-S",ylab="Est. Alpha-S");abline(0,1,col="blue")
plot(uvsd@estN[uvsd@betaN]~sim@betaN,xlab="True
Beta-N",ylab="Est. Beta-N");abline(0,1,col="blue")
plot(uvsd@estS[uvsd@betaS]~sim@betaS,xlab="True
Beta-S",ylab="Est. Beta-S");abline(0,1,col="blue")

#Sigma^2 effects
#Note that Sigma^2 is biased high with
#few participants and items.  This bias
#goes away with larger sample sizes.
par(mfrow=c(2,2),pch=19,pty='s')
matplot(sqrt(exp(uvsd@blockS2[,uvsd@muS])),t='l',xlab="Iteration",ylab="Mu-Sigma2")
abline(h=sqrt(exp(sim@muS2)),col="blue")
plot(uvsd@blockS2[,uvsd@thetaS],t='l')

plot(uvsd@estS2[uvsd@alphaS]~sim@alphaS2,xlab="True
Alpha-Sigma2",ylab="Est. Alpha-Sigma2");abline(0,1,col="blue")
plot(uvsd@estS2[uvsd@betaS]~sim@betaS2,xlab="True
Beta-Sigma2",ylab="Est. Beta-Sigma2");abline(0,1,col="blue")

#Look at some criteria
par(mfrow=c(2,2))
for(i in 1:4)
matplot(t(uvsd@s.crit[i,,]),t='l')
}

