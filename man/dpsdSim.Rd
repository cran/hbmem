\name{dpsdSim}
\alias{dpsdSim}
\title{Function dpsdSim}
\description{Simulates data from a hierarchical DPSD model.}
\usage{dpsdSim(I = 30, J = 200, K = 6, muN = -0.7, s2aN = 0.2, s2bN = 0.2,muS = 0, s2aS = 0.2, s2bS = 0.2, muR = qnorm(0.25), s2aR = 0.2,s2bR = 0.2, crit = matrix(rep(c(-1.6, -0.5, 0, 0.5, 1.6),each = I), ncol = (K - 1)))}

\arguments{
  \item{I}{Number of participants.}
  \item{J}{Number of items.}
  \item{K}{Number of response options.}
  \item{muN}{Mean of new-item distribution}
  \item{s2aN}{Variance of participant effects on mean of new-item distribution.}  \item{s2bN}{Variance of item effects on mean of new-item distribution.}
  \item{muS}{Mean of studied-item distribution}
  \item{s2aS}{Variance of participant effects on mean of studied-item distribution.}  \item{s2bS}{Variance of item effects on mean of studied-item distribution.}
  \item{muR}{Mean recollection, on probit space.}
  \item{s2aR}{Variance of participant effects recollection.}
  \item{s2bR}{Variance of item effects on recollection.}
  \item{crit}{Matrix of criteria (not including -Inf or Inf).  Columns
    correspond to criteria, rows correspond to participants.}
  }

\value{The function returns an internally defined "dpsdSim" structure.}
\references{See Pratte, Rouder, & Morey (2009)}
\author{Michael S. Pratte}
\seealso{hbmem}
\examples{
#Data from hiererchial model
sim=dpsdSim()
slotNames(sim)
table(sim@resp,sim@cond)

#Usefull to make data.frame for passing to functions
dat=as.data.frame(cbind(sim@subj,sim@item,sim@cond,sim@lag,sim@resp))
colnames(dat)=c("sub","item","cond","lag","resp")

table(dat$resp,dat$cond)
}
