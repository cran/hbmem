\name{samplePosNorm}
\alias{samplePosNorm}
\title{Function samplePosNorm}
\description{Samples posterior of mean parameters of the positive hierarchical
  linear normal model with a single Sigma2
$(x = N(exp(mu+alpha_i+beta_j),sigma2))$.}
\usage{samplePosNorm(sample, y, sub, item, lag, I, J, R, 
    sig2mu, a, b, met, sigma2, sampLag)}

\arguments{
\item{sample}{Block of linear model parameters from previous iteration.}
\item{y}{Vector of data}
\item{sub}{Vector of subject index, starting at zero.}
\item{item}{Vector of item index, starting at zero.}
\item{lag}{Vector of lag index, zero-centered.}
\item{I}{Number of subjects.}
\item{J}{Number of items.}
\item{R}{Total number of trials.}
\item{sig2mu}{Prior variance on the grand mean mu; usually set to some
  large number.}
\item{a}{Shape parameter of inverse gamma prior
  placed on effect variances.}
\item{b}{Rate parameter of inverse gamma prior
  placed on effect variances.  Setting both s2a AND s2b to be small
  (e.g., .01, .01) makes this an uninformative prior.}
  \item{met}{Vector of tuning parameter for metropolis-hastings
    sampling.  There is one tuning parameter for mu, each of I alphas,
    each of J betas, s2alpha,s2beta,and theta.  Those for s2alpha and
    s2beta are placeholders, as these parameters are sampled with gibbs.}
  \item{sigma2}{Variance of distribution.}
  \item{sampLag}{Logical.  Whether or not to sample the lag effect.}
}

\value{
  The function returns a list.  The first element of the list is the
  newly sampled block of parameters.  The second element contains a
  vector of 0s and 1s indicating which of the decorrelating steps were accepted.}
\references{Not Published yet}
\author{Michael S. Pratte}
\seealso{hbmem}
\examples{
library(hbmem)

I=25
J=40
R=I*J
t.sigma2=3
t.mu=.5
t.sig2alpha=.2
t.sig2beta=.6
t.alpha=rnorm(I,0,sqrt(t.sig2alpha))
t.beta =rnorm(J,0,sqrt(t.sig2beta))
t.theta=-.5
sub=rep(0:(I-1),each=J)
item=rep(0:(J-1),I)
lag=scale(rnorm(R,0,sqrt(t.sigma2)/10))

tmean=1:R
for(r in 1:R) tmean[r]=exp(t.mu+t.alpha[sub[r]+1]+t.beta[item[r]+1]+t.theta*lag[r])
y=rnorm(R,tmean,sqrt(t.sigma2))

M=1000 #Way too low for real analysis!
B=I+J+4
block=matrix(0,nrow=M,ncol=B)
met=rep(.1,B);jump=.0001
b0=rep(0,B)
keep=500:M
for(m in 2:M)
{
  tmp= samplePosNorm(block[m-1,],y,sub,item,lag,I,J,R,100,.01,.01,met,t.sigma2,1)
  block[m,]=tmp[[1]]
  b0=b0+tmp[[2]]

  if(m<keep[1])
  {
   met=met+(b0/m<.3)*-jump +(b0/m>.5)*jump
   met[met<jump]=jump
  }
}

est=colMeans(block[keep,])
b0/M

par(mfrow=c(3,2))
plot(exp(block[keep,1]),t='l',main="Mu",ylab="Mu")
abline(h=exp(t.mu),col="blue")
abline(h=mean(y),col="green")
acf(block[keep,1],main="ACF of Mu")

est.alpha=est[2:(I+1)]
plot(est.alpha~t.alpha,ylab="Est. Alpha",xlab="True Alpha");abline(0,1)
est.beta=est[(I+2):(I+J+1)]
plot(est.beta~t.beta,ylab="Est. Beta",xlab="True Beta");abline(0,1)

plot(block[keep,(I+J+4)],t='l',main="Theta",ylab="Theta")
abline(h=t.theta,col="blue")

plot(density(block[keep,(I+J+2)]),col="red",main="Posterior of Variances",xlim=c(0,1))
abline(v=t.sig2alpha,col="red")
lines(density(block[keep,(I+J+3)]),col="blue")
abline(v=t.sig2beta,col="blue")
}