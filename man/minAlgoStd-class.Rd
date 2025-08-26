\name{minAlgoStd-class}
\docType{class}
\alias{minAlgoStd-class}

\title{Class \code{"minAlgoStd"}}
\description{
Class for standard algorithm to minimize multivariate functions. By
standard, we mean algorithms with the main function and its gradient
input separately. Specifically, it follows \code{\link{optim}}, not
\code{\link{nlm}}. 
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("minAlgoStd",
  ...)}. It is generated my \code{\link{algoObj}}.
}
\section{Slots}{
  \describe{
    \item{\code{algo}:}{Object of class \code{"character"} ~~ }
    \item{\code{start}:}{Object of class \code{"character"} ~~ }
    \item{\code{fct}:}{Object of class \code{"character"} ~~ }
    \item{\code{grad}:}{Object of class \code{"character"} ~~ }
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
    \item{minFit}{\code{signature(object = "minAlgoStd")}: ... }
  }
}

\examples{
showClass("minAlgoStd")
}
\keyword{classes}
