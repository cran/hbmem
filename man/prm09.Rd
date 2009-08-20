\name{prm09}
\alias{prm09}
\docType{data}
\title{PRM09 Data}
\description{Confidence ratings data from Pratte, Rouder, \& Morey
  (2009).
}
\usage{data(prm09)}
\format{
  A flat-field data frame (each row is a trial) with the following variables
  \describe{
    \item{\code{cond}}{0=new; 1=studied}
    \item{\code{sub}}{index of subject starting at 0}
    \item{\code{item}}{index of item starting at 0}
    \item{\code{lag}}{index of lag, zero-centered}
    \item{\code{resp}}{which response was made; 0="sure new"}    
  }
}
\details{Participants studied a list of 240 words, and were then tested
  on the 240 studied and on 240 new words.  At test, participants made
  one of six confidence ratings ranging from "sure new" to "sure
  studied".  Note that to apply the models to these data the "Scond"
  variable should be set to ``cond", and the "cond" variable should
  be all zeros.  This is a backwards-compatibility issue.}

\source{Pratte, Rouder, \& Morey (2009).  Separating Mnemonic Process
  from Participant and Item Effects in the Assessment of ROC
  Asymmetries.  Journal of Experimental Psychology: Learning, Memory,
  and Cognition.}

\examples{
library(hbmem)
data(prm09)
table(prm09$resp,prm09$cond)
#Turn it into data suitable for
#analysis with HBMEM functions:
newdat=prm09
newdat$Scond=newdat$cond
newdat$cond=0
summary(newdat)
}
\keyword{datasets}
