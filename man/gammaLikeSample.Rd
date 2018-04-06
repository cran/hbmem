\name{gammaLikeSample}
\alias{gammaLikeSample}
\title{Function gammaLikeSample}
\description{Runs MCMC for the hierarchical Gamma Likelihood model}
\usage{gammaLikeSample(dat, M = 10000, keep = (M/10):M, getDIC = TRUE,
shape=2,jump=.005)}

\arguments{
  \item{dat}{Data frame that must include variables
    cond,sub,item,lag,resp.  Indexes for cond, sub, item, and respone
    must start at zero and have no gapes (i.e., no skipped subject
    numbers). Lags must be zero-centered.}
  \item{M}{Number of MCMC iterations.}
  \item{keep}{Which MCMC iterations should be included in estimates and
    returned.  Use keep to both get ride of burn-in, and thin chains if
    necessary}
  \item{getDIC}{Logical. should the function compute DIC value?  This
    takes a while if M is large.}
  \item{shape}{Fixed shape across both new and studied distributuions.}
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
  \item{estS2}{Not used for gamma model.}
  \item{estCrit}{Posterior means of criteria}
  \item{blockN}{Each iteration for each parameter in the new-item mean
  block.  Rows index iteration, columns index parameter.}
  \item{blockS}{Same as blockN, but for the studied-item means}
  \item{blockS2}{Not used for gamma model.}
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
  \item{b0}{Metropolis-Hastings acceptance rates for new-item
    distribution parameters.  These should be between .2 and .6.  If
    they are not, the M, keep, or jump need to be adjusted.}
   \item{b0S2}{Metropolis-Hastings acceptance rates for studied-item
    distribution parameters. }
\item{b0Crit}{Metropolis-Hastings acceptance rates for criteria. }
}

\author{Michael S. Pratte}
\seealso{hbmem}
\examples{
#This function is broken, so
#no example that works.
#make data from gamma model
if(1==0)
{
library(hbmem)
sim=gammaLikeSim(I=50,J=400,muS=log(.5),s2aS=0,s2bS=0)
dat=as.data.frame(cbind(sim@subj,sim@item,sim@cond,sim@Scond,sim@lag,sim@resp))
colnames(dat)=c("sub","item","cond","Scond","lag","resp")
dat$lag=0

table(dat$resp,dat$Scond)
M=5000
keep=500:M
gamma=gammaLikeSample(dat,M=M,keep=keep,jump=.001)

par(mfrow=c(2,3),pch=19,pty='s')
matplot(exp(gamma@blockS[,gamma@muS]),t='l',xlab="Iteration",ylab="Mu-S")
abline(h=exp(sim@muS),col="blue")
#Estimates of Alpha as function of true values
plot(gamma@estS[gamma@alphaS]~sim@alphaS,xlab="True
Alpha-S",ylab="Est. Alpha-S");abline(0,1,col="blue")
#Estimates of Beta as function of true values
plot(gamma@estS[gamma@betaS]~sim@betaS,xlab="True
Beta-S",ylab="Est. Beta-S");abline(0,1,col="blue")

#Look at some criteria
for(i in 1:3){
matplot(t(exp(gamma@s.crit[i,2:7,])),t='l')
abline(h=sim@crit[i,])
}

gamma@estS[c(gamma@s2alphaS,gamma@s2betaS)]
}

}








