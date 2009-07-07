\name{dpsdProbs}
\alias{dpsdProbs}
\title{Function dpsdProbs}
\description{Returns the probability of making confidence ratings given
  parameters of DPSD.}
\usage{dpsdProbs(r,d,crit)}

\arguments{
  \item{r}{Probability of recollection.}
  \item{d}{Mean of the signal-detection distribution.  In the common
    parameterization of the model, this would be zero for new-item
    trials, and d' for studied-item trials.  In the PRM09
    parameterization, these are dn and ds for new and studied-item
    trials, respectively.}
  \item{crit}{Criteria (not including -Inf or Inf).}
 }
  \details{For new-item trials, simply set r=0.}
\value{The function returns the probability of making each response for
  the paramters given.}
\references{See Pratte, Rouder, & Morey (2009)}
\author{Michael S. Pratte}
\seealso{hbmem}
\examples{
#Low r
dpsdProbs(.2,1,c(-1,-.5,0,.5,1))  #studied
dpsdProbs(0,-1,c(-1,-.5,0,.5,1)) #new

#High r
dpsdProbs(.6,1,c(-1,-.5,0,.5,1))  #studied
dpsdProbs(0,-1,c(-1,-.5,0,.5,1)) #new
}
