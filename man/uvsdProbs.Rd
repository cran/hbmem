\name{uvsdProbs}
\alias{uvsdProbs}
\title{Function uvsdProbs}
\description{Returns the probability of making confidence ratings given
  parameters of UVSD.}
\usage{uvsdProbs(mean,sd,bounds)}

\arguments{
  \item{mean}{Mean of the signal-detection distribution.  In the common
    parameterization of the model, this would be zero for new-item
    trials, and d' for studied-item trials.  In the PRM09
    parameterization, these are dn and ds for new and studied-item
    trials, respectively.}
  \item{sd}{Standard deviation of the distribution}
  \item{bounds}{Criteria (not including -Inf or Inf).}
 }
\value{The function returns the probability of making each response for
  the paramters given.}
\references{See Pratte, Rouder, & Morey (2009)}
\author{Michael S. Pratte}
\seealso{hbmem}
\examples{
uvsdProbs(-1,1,c(-1,-.5,0,.5,1)) #New condition
uvsdProbs(1,1.3,c(-1,-.5,0,.5,1)) #Studied condition
}
