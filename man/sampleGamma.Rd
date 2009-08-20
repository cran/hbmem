\name{sampleGamma}
\alias{sampleGamma}
\title{Function sampleGamma}
\description{Samples posterior of mean parameters of the hierarchical
  linear model on the log scale parameter of a gamma distributuion. Usually 
used within an MCMC loop.} 

\usage{sampleGamma(sample, y, cond,subj, item,
lag,N,I,J,R,ncond,nsub,nitem,s2mu, s2a, s2b, met, shape,
sampLag,pos=FALSE)}

\arguments{
\item{sample}{Block of linear model parameters from previous iteration.}
\item{y}{Vector of data}
\item{cond}{Vector fo condition index,starting at zero.}
\item{subj}{Vector of subject index, starting at zero.}
\item{item}{Vector of item index, starting at zero.}
\item{lag}{Vector of lag index, zero-centered.}
\item{N}{Numer of conditions.}
\item{I}{Number of subjects.}
\item{J}{Number of items.}
\item{R}{Total number of trials.}
\item{ncond}{Vector of length (N) containing number of trials per 
condition.}
\item{nsub}{Vector of length (I) containing number of trials per each
  subject.}
\item{nitem}{Vector of length (J) containing number of trials per each
  item.}
\item{s2mu}{Prior variance on the grand mean mu; usually set to some
  large number.}
\item{s2a}{Shape parameter of inverse gamma prior placed on effect 
variances.}
\item{s2b}{Rate parameter of inverse gamma prior
  placed on effect variances.  Setting both s2a AND s2b to be small
  (e.g., .01, .01) makes this an uninformative prior.}
  \item{met}{Vector of tuning parameter for metropolis-hastings 
steps.  Here, all sampling (except variances of alpha and beta) and 
decorrelating steps utilize the M-H sampling algorithm.  This hould 
be adjusted so that .2 < b0 < .6.}
  \item{shape}{Single shape of Gamma distribution.}
  \item{sampLag}{Logical.  Whether or not to sample the lag effect.}
\item{pos}{Logical.  If true, the model on scale is 1+exp(mu + alpha + beta).  That is, the scale is always greater than one. }
}

\value{
  The function returns a list.  The first element of the list is the
  newly sampled block of parameters.  The second element contains a
  vector of 0s and 1s indicating which of the decorrelating steps were accepted.}

\author{Michael S. Pratte}
\seealso{hbmem}
\examples{
library(hbmem)
N=2
shape=2
I=30
J=50
R=I*J
#make some data
mu=log(c(1,2))
alpha=rnorm(I,0,.2)
beta=rnorm(J,0,.2)
theta=-.001
cond=sample(0:(N-1),R,replace=TRUE)
subj=rep(0:(I-1),each=J)
item=NULL
for(i in 1:I)
item=c(item,sample(0:(J-1),J,replace=FALSE))
lag=rnorm(R,0,100)
lag=lag-mean(lag)
resp=1:R
for(r in 1:R)
{
  scale=1+exp(mu[cond[r]+1]+alpha[subj[r]+1]+beta[item[r]+1]+theta*lag[r])
  resp[r]=rgamma(1,shape=shape,scale=scale)
}

ncond=table(cond)
nsub=table(subj)
nitem=table(item)

M=500
keep=200:M
B=N+I+J+3
s.block=matrix(0,nrow=M,ncol=B)
met=rep(.08,B)
b0=rep(0,B)
jump=.0005
for(m in 2:M)
{
tmp=sampleGamma(s.block[m-1,],resp,cond,subj,item,lag,
N,I,J,R,ncond,nsub,nitem,5,.01,.01,met,2,1,pos=TRUE)
s.block[m,]=tmp[[1]]
b0=b0 + tmp[[2]]
#Auto-tuning of metropolis decorrelating steps 
if(m>20 & m<min(keep))
  {
    met=met+(b0/m<.4)*rep(-jump,B) +(b0/m>.6)*rep(jump,B)
    met[met<jump]=jump 
  }
if(m==min(keep)) b0=rep(0,B)
}

b0/length(keep) #check acceptance rate

hbest=colMeans(s.block[keep,])

par(mfrow=c(2,2),pch=19,pty='s')
matplot(s.block[keep,1:N],t='l')
abline(h=mu,col="green")
acf(s.block[keep,1])
plot(hbest[(N+1):(I+N)]~alpha)
abline(0,1,col="green")
plot(hbest[(I+N+1):(I+J+N)]~beta)
abline(0,1,col="green")



#variance of participant effect
mean(s.block[keep,(N+I+J+1)])
#variance of item effect
mean(s.block[keep,(N+I+J+2)])
#estimate of lag effect
mean(s.block[keep,(N+I+J+3)])

}


