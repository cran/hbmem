\keyword{internal}
\name{getPred}
\alias{getPred}
\title{Function getPred}
\description{getPred is a helper function that takes a block and its
  indexes, and returns the predicted value from that block for each
  trial.  It is the same thing as taking $(XtX)^(-1) Xty$  where X is the
  design matrix for 1 (mu),sub, item, and lag; y is the vector of effect
  estimates, block.}
\usage{getPred(block,cond,sub,item,lag,N,I,J,R)}

\arguments{
\item{block}{Block of parameter values starting with the grand mean,
  followed by participant effects, and finially item effects.}
\item{cond}{Vector of condition index, starting at zero.}
\item{sub}{Vector of subject index, starting at zero.}
\item{item}{Vector of item index, starting at zero.}
\item{lag}{Vector of lag index, zero-centered.}
\item{N}{Number of conditions.}
\item{I}{Number of subjects.}
\item{J}{Number of items.}
\item{R}{Total number of trials.}
}

\value{Returns predicted means for each trial.  }
\author{Michael S. Pratte}
\seealso{hbmem}
\examples{
library(hbmem)
sim=uvsdSim()
N=length(levels(as.factor(sim@cond[sim@Scond==0])))
I=length(levels(as.factor(sim@subj)))
J=length(levels(as.factor(sim@item)))
RN=sum(sim@Scond==0)

block=c(sim@muN,sim@alphaN,sim@betaN,0,0,-.001)

means=getPred(block,sim@cond[sim@Scond==0],sim@subj[sim@Scond==0],
sim@item[sim@Scond==0],sim@lag[sim@Scond==0],N,I,J,RN)

multi=sim@resp[sim@Scond==0] 

plot(multi~means,ylab="Response to New Item",xlab="Mean of New-Item Distribution")


}
