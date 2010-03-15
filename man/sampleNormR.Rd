\name{sampleNormR}
\alias{sampleNormR}
\title{Function sampleNormR}
\description{Samples posterior of mean parameters of the hierarchical
  linear normal model with the effects a linear function of some
  other variable.}
\usage{sampleNormR(sample, phi,blockD,y,subj, item, lag, I, J, R, nsub, nitem,
    s2mu, s2a, s2b, meta, metb, sigma2, sampLag)}

\arguments{
\item{sample}{Block of linear model parameters from previous iteration.}
\item{y}{Vector of data}
\item{phi}{Vector of linear slopes on effects.}
\item{blockD}{Block of parameters that will serve as the means of random effects}
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
  newly sampled block of parameters.  The THIRD element contains a
  vector of 0s and 1s indicating which of the decorrelating steps were accepted.}
\references{Not published yet.}
\author{Michael S. Pratte}
\seealso{hbmem}
\examples{
library(hbmem)

I=50
J=100
M=10
B=I+J+4
mu=.5
muS2=0
s2a=.2
s2b=.2
s2aS2=0
s2bS2=0

phi=c(.2,.08)
blockD=rep(0,B)
blockD[2:(I+1)]=rnorm(I,0,.5)
blockD[(I+2):(I+J+1)]=rnorm(J,0,.5)

    R = I * J
    alpha = rnorm(I, phi[1]*blockD[2:(I+1)], sqrt(s2a))
    beta =  rnorm(J, phi[2]*blockD[(I+2):(I+J+1)], sqrt(s2b))
    alphaS2 = rnorm(I, 0, sqrt(s2aS2))
    betaS2 = rnorm(J, 0, sqrt(s2bS2))
    subj = rep(0:(I - 1), each = J)
    item = rep(0:(J - 1), I)
    lag = rep(0, R)
    resp = 1:R
    for (r in 1:R) {
        mean = mu + alpha[subj[r] + 1] + beta[item[r] + 1]
        sd = sqrt(exp(muS2 + alphaS2[subj[r] + 1] + betaS2[item[r] + 1]))
        resp[r] = rnorm(1, mean, sd)
    }
    sim=(as.data.frame(cbind(subj, item, lag, resp)))



blockR=matrix(0,M,B)
blockR[1,c(I+J+2,I+J+3)]=c(.1,.1)
met=c(.1,.1)
b0=c(0,0)

for(m in 2:M)
  {
tmp=sampleNormR(blockR[m-1,],phi,blockD,sim$resp,sim$subj,sim$item,sim$lag,I,J,I*J,table(sim$sub),table(sim$item),10,.01,.01,met[1],met[2],1,1)
blockR[m,]=tmp[[1]]
b0=b0+tmp[[3]]
}

est=colMeans(blockR)

par(defpar(2,3))
plot(blockR[,1],t='l')
abline(h=mu,col="blue")
plot(blockR[,I+J+2],t='l')
abline(h=s2a,col="blue")
plot(blockR[,I+J+3],t='l')
abline(h=s2b,col="blue")

plot(est[2:(I+1)]~alpha);abline(0,1,col="blue")
plot(est[(I+2):(I+J+1)]~beta);abline(0,1,col="blue")

#Compare estimates from regular normal ones:

s.block=matrix(0,nrow=M,ncol=B)
met=c(.1,.1);b0=c(0,0)
for(m in 2:M)
{
tmp=sampleNorm(s.block[m-1,],sim$resp,rep(0,length(sim$resp)),sim$subj,sim$item,sim$lag,1,I,J,R,R,table(sim$subj),
table(sim$item),100,.01,.01,met[1],met[2],1,1)
s.block[m,]=tmp[[1]]
b0=b0 + tmp[[2]]
}

est2=colMeans(s.block)

par(defpar(1,2))
plot(est[2:(I+1)]~est2[2:(I+1)]);abline(0,1,col="blue")
plot(est[(I+2):(I+J+1)]~est2[(I+2):(I+J+1)]);abline(0,1,col="blue")



}