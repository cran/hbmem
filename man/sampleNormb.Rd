\name{sampleNormb}
\alias{sampleNormb}
\title{Function sampleNormb}
\description{Same as sampleNorm, but assumes an additive model on sigma2,
  and takes the block of sigma2 parameters as argument}
\usage{sampleNormb(sample,y,cond,subj,item,lag,N,I,J,R,ncond,nsub,nitem,
s2mu,s2a,s2b,meta,metb,blockSigma2,sampLag=1,Hier=1)}

\arguments{
\item{sample}{Block of linear model parameters from previous iteration.}
\item{y}{Vector of data}
\item{cond}{Vector of condition index, starting at zero.}
\item{subj}{Vector of subject index, starting at zero.}
\item{item}{Vector of item index, starting at zero.}
\item{lag}{Vector of lag index, zero-centered.}
\item{N}{Number of conditions.}
\item{I}{Number of subjects.}
\item{J}{Number of items.}
\item{R}{Total number of trials.}
\item{ncond}{Vector of length (N) containing number of trials per each
  condition.}
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
  \item{blockSigma2}{Block of parameters for Sigma2 (on log scale).
    Like all blocks, first element is the overall mean, followed by
    participant effects and then item effects.}
  \item{sampLag}{Logical.  Whether or not to sample the lag effect.}
\item{Hier}{Locial.  If TRUE then effect variances are estimated from
  data.  If false, then these values are fixed to whatever is in the
  s2alpha and s2beta slots of sample.  This value should always be
  TRUE unless you know what you are doing.}
}

\value{
  The function returns a list.  The first element of the list is the
  newly sampled block of parameters.  The second element contains a
  vector of 0s or 1s indicating which of the decorrelating steps were accepted.}
\author{Michael S. Pratte}
\seealso{hbmem,sampleSig2b}
\examples{
library(hbmem)
N=2
I=50
J=200
B=N+I+J+3
R = I * J

mu=c(3,5)
muS2=log(c(1,2))
alpha = rnorm(I, 0, sqrt(.2))
beta = rnorm(J, 0, sqrt(.2))
alphaS2 = rnorm(I, 0, sqrt(.2))
betaS2 = rnorm(J, 0, sqrt(.2))
cond=sample(0:(N-1),R,replace=TRUE)
subj = rep(0:(I - 1), each = J)
item = rep(0:(J - 1), I)
lag = rep(0, R)
lag=runif(R,-500,500)
lag=lag-mean(lag)
resp = 1:R
for (r in 1:R) {
    mean = mu[cond[r] + 1] + alpha[subj[r] + 1] + beta[item[r] + 1]
    sd = sqrt(exp(muS2[cond[r]+1] + alphaS2[subj[r] + 1] +
betaS2[item[r] + 1] + .005*lag[r]))
    resp[r] = rnorm(1, mean, sd)
}
sim=(as.data.frame(cbind(cond,subj, item, lag, resp)))
attach(sim)
plot(resp~lag)

########MCMC SETUP##########
blockS=blockS2=matrix(0,nrow=10,ncol=B)
blockS[,B-1]=blockS[,B-2]=blockS2[,B-1]=blockS2[,B-2]=.5
b0mean=c(0,0)
b0S2=rep(0,B)
met=rep(.01,B)
jump=.0001
ncond=table(cond)
nsub=table(subj)
nitem=table(item)

for(m in 2:10) #way to low for real analysis
  {
    tmp=sampleNormb(blockS[m-1,],resp,cond,subj,item,lag,N,I,J,I*J,
ncond,nsub,nitem,10,.01,.01,.02,.005,blockS2[m-1,],1,1)
    blockS[m,]=tmp[[1]]
    b0mean=b0mean+tmp[[2]]
    
    tmp=sampleSig2b(blockS2[m-1,],resp,cond,subj,item,lag,N,I,J,I*J,
ncond,nsub,nitem,10,.01,.01,met,blockS[m,],1,1)
    blockS2[m,]=tmp[[1]]
    b0S2=b0S2+tmp[[2]]
   if(m<10) met=met+(b0S2/m<.3)*-jump +(b0S2/m>.5)*jump
    met[met<jump]=jump
#met[B]=.0001
  }
b0mean/m
b0S2/m

s=colMeans(blockS)
s2=colMeans(blockS2)

par(mfrow=c(3,3))
matplot(blockS[,1:N],t='l')
abline(h=mu)
plot(s[(N+1):(I+N)]~alpha);abline(0,1)
plot(s[(I+N+1):(I+J+N)]~beta);abline(0,1)

matplot(blockS2[,1:N],t='l')
abline(h=muS2)
plot(s2[(N+1):(I+N)]~alphaS2);abline(0,1)
plot(s2[(I+N+1):(I+N+J)]~betaS2);abline(0,1)

plot(blockS2[,B-2],t='l')
plot(blockS2[,B-1],t='l')
plot(blockS2[,B],t='l')

}
