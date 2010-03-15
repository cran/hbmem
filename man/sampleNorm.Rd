\name{sampleNorm}
\alias{sampleNorm}
\title{Function sampleNorm}
\description{Samples posterior of mean parameters of the hierarchical
  linear normal model with a single Sigma2.  Usually used within an MCMC
  loop.}
\usage{sampleNorm(sample, y, cond, subj, item, lag, N, I, J, R, ncond, nsub,
nitem, s2mu, s2a, s2b, meta, metb, sigma2, sampLag=TRUE,Hier=TRUE)}

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
  \item{sigma2}{Variance of distribution.}
  \item{sampLag}{Logical.  Whether or not to sample the lag effect.}
  \item{Hier}{Logical. If TRUE then effect variances are estimated
    from data.  If FALSE then these values are set to whatever
    value is in the s2alpha and s2beta slots of sample.  This
    should always be set to TRUE.}
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
N=2
t.mu=c(1,2)
I=20
J=50
R=I*J
#make some data
tmp=normalSim(N=N,I=I,J=J,mu=t.mu,s2a=2,s2b=2,muS2=log(1),s2aS2=0,s2bS2=0)
dat=tmp[[1]]
t.alpha=tmp[[2]]
t.beta=tmp[[3]]

ncond=table(dat$cond)
nsub=table(dat$sub)
nitem=table(dat$item)

M=10
keep=2:M
B=N+I+J+3
s.block=matrix(0,nrow=M,ncol=B)
met=c(.1,.1);b0=c(0,0)
jump=.001
for(m in 2:M)
{
tmp=sampleNorm(s.block[m-1,],dat$resp,dat$cond,dat$subj,dat$item,dat$lag,
N,I,J,R,ncond,nsub,nitem,5,.01,.01,met[1],met[2],1,1,1)
s.block[m,]=tmp[[1]]
b0=b0 + tmp[[2]]


#Auto-tuning of metropolis decorrelating steps 
if(m>20 & m<min(keep))
  {
    met=met+(b0/m<.2)*c(-jump,-jump) +(b0/m>.3)*c(jump,jump)
    met[met<jump]=jump
  }
}

b0/M #check acceptance rate

hbest=colMeans(s.block[keep,])

par(mfrow=c(2,2),pch=19,pty='s')
matplot(s.block[keep,1:N],t='l')
abline(h=t.mu,col="green")
abline(h=tapply(dat$resp,dat$cond,mean),col="orange")
acf(s.block[keep,1])
plot(hbest[(N+1):(I+N)]~t.alpha)
abline(0,1,col="green")
plot(hbest[(I+N+1):(I+J+N)]~t.beta)
abline(0,1,col="green")



#variance of participant effect
mean(s.block[keep,(N+I+J+1)])
#variance of item effect
mean(s.block[keep,(N+I+J+2)])
#estimate of lag effect
mean(s.block[keep,(N+I+J+3)])

}


