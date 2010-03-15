\name{dpsdSim}
\alias{dpsdSim}
\title{Function dpsdSim}
\description{Simulates data from a hierarchical DPSD model.}
\usage{dpsdSim(NN=2,NS=1,I=30,J=200,K=6,muN=c(-.7,-.5),s2aN=.2,s2bN=.2,
muS=0,s2aS=.2,s2bS=.2,muR=qnorm(.25),s2aR=.2,s2bR=.2,
crit=matrix(rep(c(-1.6,-.5,0,.5,1.6),each=I),ncol=(K-1)))}

\arguments{
  \item{NN}{Number of new-item conditions.}
  \item{NS}{Number of studied-item conditions.}
  \item{I}{Number of participants.}
  \item{J}{Number of items.}
  \item{K}{Number of response options.}
  \item{muN}{Mean of new-item distribution.  If there are more than one new-item conditions this is a vector of means with length equal to NN.}
  \item{s2aN}{Variance of participant effects on mean of new-item distribution.}  
  \item{s2bN}{Variance of item effects on mean of new-item distribution.}
  \item{muS}{Mean of studied-item distribution.  If there are more than
    new-item conditions this is a vector of means with length equal to NNone studied-item conditions this is a vector of means with length equal to NS.}
  \item{s2aS}{Variance of participant effects on mean of studied-item distribution.}  
  \item{s2bS}{Variance of item effects on mean of studied-item distribution.}
  \item{muR}{Mean recollection, on probit space.}
  \item{s2aR}{Variance of participant effects recollection.}
  \item{s2bR}{Variance of item effects on recollection.}
  \item{crit}{Matrix of criteria (not including -Inf or Inf).  Columns
    correspond to criteria, rows correspond to participants.}
  }

\value{The function returns an internally defined "dpsdSim" structure.}
\references{See Pratte, Rouder, & Morey (2009)}
\author{Michael S. Pratte}
\keyword{models}
\seealso{hbmem}
\examples{
library(hbmem)
#Data from hiererchial model
sim=dpsdSim()
slotNames(sim)
#Scond indicates studied/new
#cond indicates which condition (e.g., deep/shallow)

table(sim@resp,sim@Scond,sim@cond)

#Usefull to make data.frame for passing to functions
dat=as.data.frame(cbind(sim@subj,sim@item,sim@Scond,sim@cond,sim@lag,sim@resp))
colnames(dat)=c("sub","item","Scond","cond","lag","resp")

table(dat$resp,dat$Scond,dat$cond)
}
