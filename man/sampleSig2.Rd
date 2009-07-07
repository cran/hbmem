\name{sampleSig2}
\alias{sampleSig2}
\title{Function sampleSig2}
\description{Samples posterior of mean of the variance of a normal
  distibution which has an additive structure on the mean, and a single
  variance for all values.  Usually used within MCMC loop.}
\usage{sampleSig2(sig2,block,y,sub,item,lag,I, J, R,a,b)}

\arguments{
\item{sig2}{Sample of sig2 from previous iteration.}
\item{block}{Vector of parameters for mean of distribution}
\item{y}{Vector of data}
\item{sub}{Vector of subject index, starting at zero.}
\item{item}{Vector of item index, starting at zero.}
\item{lag}{Vector of lag index, zero-centered.}
\item{I}{Number of subjects.}
\item{J}{Number of items.}
\item{R}{Total number of trials.}
\item{a}{Shape parameter for inverse gamma prior on Sigma2.}
\item{b}{Rate parameter for inverse gamma prior on Sigma2.  Setting 'a'
  and 'b' to small values (e.g., .01, .01) makes the prior non-informative.}
}

\value{The function returns a new sample of Sigma2.}
\author{Michael S. Pratte}
\seealso{hbmem}
\examples{
library(hbmem)
true.mean=0
true.sigma2=1.4
I=1
J=1
R=100
sub=rep(0,R)
item=rep(0,R)
lag=rep(0,R)

#make some data
dat=rnorm(R,true.mean,sqrt(true.sigma2))
true.block=c(true.mean,rep(0,(I+J+3)))
a=b=.01

M=1000
s.sigma2=rep(1,M)

for(m in 2:M)
{
s.sigma2[m]=sampleSig2(s.sigma2[m-1],true.block,dat,sub,item,lag,I, J, R,a,b)
}

par(mfrow=c(1,3),pty='s')
hist(s.sigma2,prob=TRUE)
lines(density(s.sigma2),col="red")

plot(s.sigma2,t='l')
abline(h=true.sigma2,col="blue")
abline(h=mean(s.sigma2),col="green") #post mean

acf(s.sigma2)

#95 credible interval
quantile(s.sigma2,c(.025,.975))
}
