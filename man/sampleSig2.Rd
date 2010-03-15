\name{sampleSig2}
\alias{sampleSig2}
\title{Function sampleSig2}
\description{Samples posterior of the variance of a normal
  distibution which has an additive structure on the mean, and a single
  variance for all values.  Usually used within MCMC loop.}
\usage{sampleSig2(sig2,block,y,cond,sub,item,lag,N,ncond,I,J,a,b)}

\arguments{
\item{sig2}{Sample of sig2 from previous iteration.}
\item{block}{Vector of parameters for mean of distribution}
\item{y}{Vector of data}
\item{cond}{Vector that indexs condition (e.g., deep vs. shallow)}
\item{sub}{Vector of subject index, starting at zero.}
\item{item}{Vector of item index, starting at zero.}
\item{lag}{Vector of lag index, zero-centered.}
\item{N}{Number of conditions.}
\item{ncond}{Number of trials per condition.}
\item{I}{Number of subjects.}
\item{J}{Number of items.}
\item{a}{Shape parameter for inverse gamma prior on Sigma2.}
\item{b}{Rate parameter for inverse gamma prior on Sigma2.  Setting 'a'
  and 'b' to small values (e.g., .01, .01) makes the prior non-informative.}
}

\value{The function returns a new sample of Sigma2.}
\author{Michael S. Pratte}
\seealso{hbmem}
\examples{
library(hbmem)
true.mean=c(0,0)
true.sigma2=c(10,20)
N=2
I=1
J=1
R=1000
cond=rep(0:1,R/2)
ncond=table(cond)
sub=rep(0,R)
item=rep(0,R)
lag=rep(0,R)

#make some data
dat=rnorm(R,true.mean[cond+1],sqrt(true.sigma2[cond+1]))
true.block=c(true.mean,rep(0,(I+J+3)))

a=b=.01

M=10
s.sigma2=matrix(1,M,N)

for(m in 2:M)
{
s.sigma2[m,]=sampleSig2(s.sigma2[m-1,],true.block,dat,cond,sub,item,lag,N,ncond,I, J,a,b)
}

par(mfrow=c(1,1),pty='s')

matplot(s.sigma2,t='l')

abline(h=true.sigma2,col="blue")
abline(h=colMeans(s.sigma2),col="green") #post mean

}
