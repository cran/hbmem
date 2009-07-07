\name{rtnorm}
\alias{rtnorm}
\title{Function rtnorm}
\description{Returns random samples from a truncated normal distribution.}
\usage{rtnorm(N,mu,sigma,a,b)}

\arguments{
  \item{N}{Number of samples to return.}
  \item{mu}{A vector of length N that contains distribution means for
    each draw.}
  \item{sigma}{A vector of length N that contains distribution standard
deviations for each draw.}
  \item{a}{Vector of length N of lower truncation points; may be -Inf.}
  \item{b}{Vector of length N of upper truncation point; may be Inf.} 
}
  \details{This function is currently unstable for drawing from regions
    with extremely low probabilities.  If this happens is should print
    a warning, and return a draw from a uniform distribution between a
    and b.  See example below for how to break it.}
\value{Returns 'N' random draws.}
\author{Michael S. Pratte}
\seealso{hbmem}
\examples{
library(hbmem)
#Draw one
rtnorm(1,0,1,0,.2)

#Draw 50
N=500
mu=rep(0,N)
sigma=rep(1,N)
a=rep(1,N)
b=rep(2,N)
x=rtnorm(N,mu,sigma,a,b)
hist(x)

#Break it
rtnorm(1,0,1,1000,1001)
}
