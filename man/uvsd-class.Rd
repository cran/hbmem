\name{uvsd-class}
\docType{class}
\alias{uvsd-class}

\title{Class "uvsd" ~~~ }
\description{Object that results from running hierarchical UVSD analysis.}

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
    \item{\code{estS2}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{estCrit}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{blockN}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{blockS}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{blockS2}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{s.crit}:}{Object of class \code{"array"} ~~ }
    \item{\code{pD}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{DIC}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{M}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{keep}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{b0}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{b0S2}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{b0Crit}:}{Object of class \code{"numeric"} ~~ }
  }
}
\section{Methods}{
No methods defined with class "uvsd" in the signature.
}
\seealso{uvsdSample()}
\examples{
showClass("uvsd")
}
\keyword{classes}
