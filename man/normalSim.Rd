\name{normalSim}
\alias{normalSim}
\title{Function normalSim}
\description{Simulates data from a hierarchical linear normal model.}
\usage{normalSim(I=30,J=300,mu=-.5,s2a=.2,s2b=.2,muS2=0,s2aS2=.2,s2bS2=.2)}

\arguments{
  \item{I}{Number of participants.}
  \item{J}{Number of items.}
  \item{mu}{Grand mean}
  \item{s2a}{Variance of subject effect on the mean}
  \item{s2b}{Variance of item effect on the mean}
  \item{muS2}{Overall variance of data on log scale}
  \item{s2aS2}{Variance of subject effect on variance}
  \item{s2bS2}{Variance of item effect on variance}
}

\value{The function returns a data frame with subject (subj), item, lag,
  and response (resp) columns.  Lag is a vector of zeros (i.e., no lag effect).}
\author{Michael S. Pratte}
\seealso{hbmem}
\examples{
library(hbmem)
I=20
J=50
R=I*J
dat=normalSim(I=I,J=J,mu=10,s2a=1,s2b=1,muS2=log(1),s2aS2=0,s2bS2=0)
summary(dat)
}
