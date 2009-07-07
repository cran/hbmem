\name{sampleNorm}
\alias{sampleNorm}
\title{Function sampleNorm}
\description{Samples posterior of mean parameters of the hierarchical
  linear normal model with a single Sigma2.  Usually used within an MCMC
loop.}
\usage{sampleNorm(sample, y, subj, item, lag, I, J, R, nsub, nitem,
    s2mu, s2a, s2b, meta, metb, sigma2, sampLag)}

\arguments{
\item{sample}{Block of linear model parameters from previous iteration.}
\item{y}{Vector of data}
\item{subj}{Vector of subject index, starting at zero.}
\item{item}{Vector of item index, starting at zero.}
\item{lag}{Vector of lag index, zero-centered.}
\item{I}{Number of subjects.}
\item{J}{Number of items.}
\item{R}{Total number of trials.}
\item{nsub}{Vector of length (I) containing number of trials per each
  subject.}
\item{nitem}{Vector of length (J) containing number of trials per each
  item.}
\item{s2mu}{Prior variance on the grand mean mu; usually set to some
  large number.}
\item{s2a}{Shape parameter of inverse gamma prior
  placed on effect variances.}
\item{s2b}{Rate parameter of inverse gamma prior
  placed on effect variances.  Setting both s2a AND s2b to be small
  (e.g., .01, .01) makes this an uninformative prior.}
  \item{meta}{Matrix of tuning parameter for metropolis-hastings decorrelating
    step on mu and alpha.  This hould be adjusted so that .2 < b0 < .6.}
  \item{metb}{Tunning parameter for decorrelating step on alpha and beta.}
  \item{sigma2}{Variance of distribution.}
  \item{sampLag}{Logical.  Whether or not to sample the lag effect.}
}

\value{
  The function returns a list.  The first element of the list is the
  newly sampled block of parameters.  The second element contains a
  vector of 0s and 1s indicating which of the decorrelating steps were accepted.}
\references{See Pratte, Rouder, & Morey (2009)}
\author{Michael S. Pratte}
\seealso{hbmem}
\examples{
library(hbmem)
I=20
J=50
R=I*J
#make some data
dat=normalSim(I=I,J=J,mu=10,s2a=1,s2b=1,muS2=log(1),s2aS2=0,s2bS2=0)
nsub=table(dat$sub)
nitem=table(dat$item)

M=2000
keep=200:M
B=I+J+4
s.block=matrix(0,nrow=M,ncol=B)
met=c(.1,.1);b0=c(0,0)
for(m in 2:M)
{
tmp=sampleNorm(s.block[m-1,],dat$resp,dat$subj,dat$item,dat$lag,I,J,R,nsub,nitem,100,.01,.01,met[1],met[2],1,1)
s.block[m,]=tmp[[1]]
b0=b0 + tmp[[2]]
}

hbest=colMeans(s.block[keep,])
estAlpha=tapply(dat$resp,dat$subj,mean) - mean(dat$resp)
estBeta=tapply(dat$resp,dat$item,mean) - mean(dat$resp)

par(mfrow=c(2,3),pch=19,pty='s')
plot(s.block[keep,1],t='l')
abline(h=mean(dat$resp),col="green")
plot(hbest[2:(I+1)]~estAlpha)
abline(0,1,col="green")
plot(hbest[(I+2):(I+J+1)]~estBeta)
abline(0,1,col="green")

#variance of participant effect
hist(s.block[keep,(I+J+2)])
#variance of item effect
hist(s.block[keep,(I+J+3)])
#estimate of lag effect 
hist(s.block[keep,(I+J+4)])
}
