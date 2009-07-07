\name{dpsd-class}
\docType{class}
\alias{dpsd-class}

\title{Class "dpsd" ~~~ }
\description{Object that is output from performing estimation of
  hierarchical DPSD model; contains all estimates, chains, etc.}
\section{Slots}{
	 \describe{
    \item{\code{mu}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{alpha}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{beta}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{s2alpha}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{s2beta}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{theta}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{estN}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{estS}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{estR}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{estCrit}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{blockN}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{blockS}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{blockR}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{s.crit}:}{Object of class \code{"array"} ~~ }
    \item{\code{pD}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{DIC}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{M}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{keep}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{b0}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{b0Crit}:}{Object of class \code{"numeric"} ~~ }
  }
}
\section{Methods}{
No methods defined with class "dpsd" in the signature.
}
\author{Michael Pratte}

\seealso{dpsdSample}
\examples{
#See help for 'hbmem'
}
\keyword{classes}
