\name{dpsdPosLogLike}
\alias{dpsdPosLogLike}
\title{function dpsdPosLogLike}
\description{Computes log likelihood for DPSD model with positive d'}
\usage{
dpsdPosLogLike(R, NN, NS, I, JN, JS, K, dat, cond, Scond, sub, item,lag, blockN, blockD, blockR, crit)
}
\arguments{
\item{R}{Total number of trials.}
\item{NN}{Number of new-item conditions.}
\item{NS}{Number of studied-item conditions.}
\item{I}{Number of subjects.}
\item{JN}{Number of items in new condition.}
\item{JS}{Number of items in studied condition.}
\item{K}{Number of response options.}
\item{dat}{Vector of responses, ranging from 0:(K-1).}
\item{cond}{Vector of condition index.}
\item{Scond}{Vector of new-studied condition index; 0=new, 1=studied.}
\item{sub}{Vector of subject index, starting at 0 with no missing
  subject numbers.}
\item{item}{Vector of item index, starting at 0 with no missing
  item numbers.}
\item{lag}{Vector of lag index.}
\item{blockN}{Block of parameters for new-item means.}
\item{blockD}{Block of parameters for dprime means.}
\item{blockR}{Block of parameters for recollection values.}
\item{crit}{VECTOR of criteria including -Inf and Inf for top and bottom
critieria, respectively.  Vector contains the (K+1) criteria for the
first subjects, followed by those for the second subject, etc.}
}
\keyword{internal}