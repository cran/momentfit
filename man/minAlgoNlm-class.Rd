\name{minAlgoNlm-class}
\docType{class}
\alias{minAlgoNlm-class}

\title{Class \code{"minAlgoNlm"}}
\description{
Class for algorithm to minimize multivariate functions that have the
same format as \code{\link{nlm}}. 
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("minAlgoNlm",
  ...)}. It is generated my \code{\link{algoObj}}.
}
\section{Slots}{
  \describe{
    \item{\code{algo}:}{Object of class \code{"character"} ~~ }
    \item{\code{start}:}{Object of class \code{"character"} ~~ }
    \item{\code{fct}:}{Object of class \code{"character"} ~~ }
    \item{\code{solution}:}{Object of class \code{"character"} ~~ }
    \item{\code{value}:}{Object of class \code{"character"} ~~ }
    \item{\code{message}:}{Object of class \code{"character"} ~~ }
    \item{\code{convergence}:}{Object of class \code{"character"} ~~ }
  }
}
\section{Extends}{
Class \code{"\linkS4class{minAlgo}"}, directly.
}
\section{Methods}{
  \describe{
    \item{minFit}{\code{signature(object = "minAlgoNlm")}: ... }
  }
}

\examples{
showClass("minAlgoNlm")
}
\keyword{classes}
