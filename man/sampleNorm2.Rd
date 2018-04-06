\name{sampleNorm2}
\alias{sampleNorm2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Function sampleNorm2}
\description{Function for sampling normal distribution with grand means
  built into participant effects (i.e., participant effects not
  zero-centered).  This is useful for sampling DPSD or UVSD models with
  no item effects and without hierarchical pooling (neither of which
  should be done!).}
\usage{
sampleNorm2(sample, y, cond, subj, item, lag, N, I, J, R,
ncond, nsub, nitem, s2mu, s2a, s2b, meta, metb, sigma2, sampLag = 1, Hier = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sample}{ ~~Describe \code{sample} here~~ }
  \item{y}{ ~~Describe \code{y} here~~ }
  \item{cond}{ ~~Describe \code{cond} here~~ }
  \item{subj}{ ~~Describe \code{subj} here~~ }
  \item{item}{ ~~Describe \code{item} here~~ }
  \item{lag}{ ~~Describe \code{lag} here~~ }
  \item{N}{ ~~Describe \code{N} here~~ }
  \item{I}{ ~~Describe \code{I} here~~ }
  \item{J}{ ~~Describe \code{J} here~~ }
  \item{R}{ ~~Describe \code{R} here~~ }
  \item{ncond}{ ~~Describe \code{ncond} here~~ }
  \item{nsub}{ ~~Describe \code{nsub} here~~ }
  \item{nitem}{ ~~Describe \code{nitem} here~~ }
  \item{s2mu}{ ~~Describe \code{s2mu} here~~ }
  \item{s2a}{ ~~Describe \code{s2a} here~~ }
  \item{s2b}{ ~~Describe \code{s2b} here~~ }
  \item{meta}{ ~~Describe \code{meta} here~~ }
  \item{metb}{ ~~Describe \code{metb} here~~ }
  \item{sigma2}{ ~~Describe \code{sigma2} here~~ }
  \item{sampLag}{ ~~Describe \code{sampLag} here~~ }
  \item{Hier}{ ~~Describe \code{Hier} here~~ }
}
\keyword{internal}