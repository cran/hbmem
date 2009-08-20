\name{uvsdLogLike}
\alias{uvsdLogLike}
\title{Function uvsdLogLike}
\description{Computes log likelihood for UVSD model}
\usage{uvsdLogLike(R,NN,NS,I,J,K,dat,cond,Scond,subj,item,lag,blockN,blockS,blockS2,crit)}

\arguments{
\item{R}{Total number of trials.}
\item{NN}{Number of new-item conditions.}
\item{NS}{Number of studied-item conditions.}
\item{I}{Number of subjects.}
\item{J}{Number of items.}
\item{K}{Number of response options.}
\item{dat}{Vector of responses, ranging from 0:(K-1).}
\item{cond}{Vector of condition index.}
\item{Scond}{Vector of new/studied condition index; 0=new, 1=studied.}
\item{subj}{Vector of subject index, starting at 0 with no missing
  subject numbers.}
\item{item}{Vector of item index, starting at 0 with no missing
  item numbers.}
\item{lag}{Vector of lag index.}
\item{blockN}{Block of parameters for new-item means.}
\item{blockS}{Block of parameters for studied-item means.}
\item{blockS2}{Block of parameters for Sigma2 values.  If there is only
  one Sigma2 for all participants and items, then the first element of
  blockS2 should contain this value, and the other elements fo blockS2 should
  be zero.}
\item{crit}{VECTOR of criteria including -Inf and Inf for top and bottom
critieria, respectively.  Vector contains the (K+1) criteria for the
first subjects, followed by those for the second subject, etc.}
}
\value{The function returns the log likelihood.}
\references{See Pratte, Rouder, & Morey (2009)}
\author{Michael S. Pratte}
\seealso{hbmem}
