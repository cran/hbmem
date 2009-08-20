\name{uvsdSim}
\alias{uvsdSim}
\title{Function uvsdSim}
\description{Simulates data from a hierarchical UVSD model.}
\usage{uvsdSim(NN = 2, NS = 1, I = 30, J = 200, K = 6, muN = c(-0.5, 
    -0.2), s2aN = 0.2, s2bN = 0.2, muS = 0.5, s2aS = 0.2, s2bS = 0.2, 
    muS2 = log(1), s2aS2 = 0, s2bS2 = 0,lagEffect = -0.001,crit = matrix(rep(c(-1.5,-0.5, 0, 0.5, 1.5), each = I), ncol = (K - 1)))}

\arguments{
  \item{NN}{Number of conditions for new words.}
  \item{NS}{Number of conditions for studied words.}
  \item{I}{Number of participants.}
  \item{J}{Number of items.}
  \item{K}{Number of response options.}
  \item{muN}{Mean of new-item distribution.  If NN is greater than 1, then muN must be a vector of length NN.}
  \item{s2aN}{Variance of participant effects on mean of new-item distribution.}  \item{s2bN}{Variance of item effects on mean of new-item distribution.}
  \item{muS}{Mean of studied-item distribution.  If NS is greater than 1, then muS must be a vector of length NS.}
  \item{s2aS}{Variance of participant effects on mean of studied-item
    distribution.}
  \item{s2bS}{Variance of item effects on mean of studied-item distribution.}
\item{lagEffect}{Magnitude of linear lag effect on both studied-item
  distribution and log(sigma2).}
  \item{muS2}{Mean variance of studied-item distribution, sigma2}
  \item{s2aS2}{Variance of participant effects sigma2.}
  \item{s2bS2}{Variance of item effects on sigma2.}
  \item{crit}{Matrix of criteria (not including -Inf or Inf).  Columns
    correspond to criteria, rows correspond to participants.}
  }

\value{The function returns an internally defined "uvsdSim" structure.}
\references{See Pratte, Rouder, & Morey (2009)}
\author{Michael S. Pratte}
\seealso{hbmem}
\examples{
library(hbmem)
#Data from hiererchial model
sim=uvsdSim() 
slotNames(sim) 
table(sim@resp,sim@Scond,sim@cond)

#Usefull to make data.frame for passing to model-fitting functions
dat=as.data.frame(cbind(sim@subj,sim@item,sim@cond,sim@Scond,sim@lag,sim@resp))
colnames(dat)=c("sub","item","cond","Scond","lag","resp")

table(dat$resp,dat$Scond,dat$cond)
}
