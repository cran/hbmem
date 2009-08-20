\name{gammaSim}
\alias{gammaSim}
\title{Function gammaSim}
\description{Simulates data from a hierarchical Gamma model.}
\usage{gammaSim(NN=1,NS=2,I=30,J=200,K=6,muN=log(.65),s2aN=.2,s2bN=.2,
muS=log(c(.8,1.2)),s2aS=.2,s2bS=.2,lagEffect=-.001,shape=2,
crit=matrix(rep(c(.3,.6,1,1.2,1.6),each=I),ncol=(K-1)))}

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
  \item{lagEffect}{Linear slope of lag effect on log of studied-item scale.}
  \item{shape}{Common shape for both new and studied distributuions.}
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
sim=gammaSim() 
slotNames(sim) 
table(sim@resp,sim@cond,sim@Scond)

#Usefull to make data.frame for passing to model-fitting functions
dat=as.data.frame(cbind(sim@subj,sim@item,sim@cond,sim@Scond,sim@lag,sim@resp))
colnames(dat)=c("sub","item","cond","Scond","lag","resp")

table(dat$resp,dat$cond,dat$Scond)
}
