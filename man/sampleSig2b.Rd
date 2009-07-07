\name{sampleSig2b}
\alias{sampleSig2b}
\title{Function sampleSig2b}
\description{Samples posterior of the variance of a normal
  distibution which has the same additive structure on the mean and the
  log of variance.  Usually used within MCMC loop, but it is broken as
  of version .1!}
\usage{sampleSig2b(sample,y,sub,item,lag,I,J,R,nsub,nitem,s2mu,s2a,s2b,met,blockMean,sampLag)}

\arguments{
\item{sample}{Previous sample of block variances.}
\item{y}{Vector of data}
\item{sub}{Vector of subject index, starting at zero.}
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
\item{s2a}{Shape parameter of inverse gamma prior placed on effect variances.}
\item{s2b}{Rate parameter of inverse gamma prior
  placed on effect variances.  Setting both s2a AND s2b to be small
  (e.g., .01, .01) makes this an uninformative prior.}
\item{met}{Vector of metropolis-hastins tuning parameters.}
\item{blockMean}{Block of parameters for the mean of the distribution.}
\item{sampLag}{Logical.  Whether or not to sample the lag effect.}
}

\value{The function returns a new sample of a block of Sigma2 paramters.}
\details{This function is for a model with an additive structure on the
  log of the variance of a normal distribuiton.  This model is under
  development, the code is buggy, and it might not even work in the
  end.}
\author{Michael S. Pratte}
\seealso{hbmem}
\examples{
#no example yet
}
