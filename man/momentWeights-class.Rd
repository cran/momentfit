\name{momentWeights-class}
\docType{class}
\alias{momentWeights-class}
\alias{momentWeights}

\title{Class \code{"momentWeights"}}
\description{
  A class to store the weighting matrix of a set of moment conditions.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("momentWeights", ...)}.
It is created my \code{\link{evalWeights}}.
}
\section{Slots}{
  \describe{
    \item{\code{w}:}{Object of class \code{"ANY"} ~~ }
    \item{\code{type}:}{Object of class \code{"character"} ~~ }
    \item{\code{wSpec}:}{Object of class \code{"list"} ~~ }
  }
}
\section{Methods}{
  \describe{
    \item{[}{\code{signature(x = "momentWeights", i = "missing", j = "missing")}: ... }
    \item{[}{\code{signature(x = "momentWeights", i = "numeric", j = "missing")}: ... }
    \item{evalGmmObj}{\code{signature(object = "momentModel", theta = "numeric", wObj = "momentWeights")}: ... }
    \item{print}{\code{signature(x = "momentWeights")}: ... }
    \item{quadra}{\code{signature(w = "momentWeights", x = "matrixORnumeric", y = "matrixORnumeric")}: ... }
    \item{quadra}{\code{signature(w = "momentWeights", x = "matrixORnumeric", y = "missing")}: ... }
    \item{quadra}{\code{signature(w = "momentWeights", x = "missing", y = "missing")}: ... }
    \item{show}{\code{signature(object = "momentWeights")}: ... }
    \item{solveGmm}{\code{signature(object = "allNLModel", wObj = "momentWeights")}: ... }
    \item{solveGmm}{\code{signature(object = "linearModel", wObj = "momentWeights")}: ... }
	 }
}

\examples{
showClass("momentWeights")
}
\keyword{classes}
