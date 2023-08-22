\name{dpsdSample}
\alias{dpsdSample}
\title{Function to fit hierarchical DPSD model to data.}
\description{Runs MCMC estimation for the hierarchical DPSD model.}
\usage{dpsdSample(dat, M = 5000, keep = (M/10):M, getDIC = TRUE,
freeCrit=TRUE,Hier=TRUE, jump=.01)}

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
  \item{freeCrit}{Logical.  If true then criteria are estimated
    separately
    for each participant.  Should be set to false if analizing only one
    participant (e.g., if averaging over subjects).} 
  \item{Hier}{Logical. If true then the variances of effects
    (e.g., item effects) are estimated from the data, i.e., effects are
    treated as random.  If false then these variances are fixed to
    2.0 (.5 for recollection effects), thus treating these effects as
    fixed.  This option is there to allow for compairson with more
    traditional approaches, and to see the effects of imposing
    hierarcical structure.  It should always be set to TRUE in real
    analysis, and is not even guaranteed to work if set to false.}
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
\keyword{models}
\seealso{hbmem}
\examples{
#In this example we generate data from EVSD, then fit it with both
#hierarchical DPSD and DPSD assuming no participant or item effects.
library(hbmem)
sim=dpsdSim(I=30,J=200)
dat=as.data.frame(cbind(sim@subj,sim@item,sim@cond,sim@Scond,sim@lag,sim@resp))
colnames(dat)=c("sub","item","cond","Scond","lag","resp")
dat$lag[dat$Scond==1]=dat$lag[dat$Scond==1]-mean(dat$lag[dat$Scond==1])

M=10 #Too low for real analysis!
keep=2:M
DPSD=dpsdSample(dat,M=M)

#Look at all parameters
par(mfrow=c(3,3),pch=19,pty='s')

matplot(DPSD@blockN[,DPSD@muN],t='l',
ylab="muN")
abline(h=sim@muN,col="blue")
plot(DPSD@estN[DPSD@alphaN]~sim@alphaN)
abline(0,1,col="blue")
plot(DPSD@estN[DPSD@betaN]~sim@betaN)
abline(0,1,col="blue")

matplot(DPSD@blockS[,DPSD@muS],t='l',
ylab="muS")
abline(h=sim@muS,col="blue")
plot(DPSD@estS[DPSD@alphaS]~sim@alphaS)
abline(0,1,col="blue")
plot(DPSD@estS[DPSD@betaS]~sim@betaS)
abline(0,1,col="blue")

matplot(pnorm(DPSD@blockR[,DPSD@muS]),t='l',
ylab="P(recollection)")
abline(h=pnorm(sim@muR),col="blue")
plot(DPSD@estR[DPSD@alphaS]~sim@alphaR)
abline(0,1,col="blue")
plot(DPSD@estR[DPSD@betaS]~sim@betaR)
abline(0,1,col="blue")
}
