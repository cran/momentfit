\name{sysMomentWeights-class}
\docType{class}
\alias{sysMomentWeights-class}

\title{Class \code{"sysMomentWeights"}}
\description{
A class to store the weighting matrix of the moment conditions from a
system of equations.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("sysMomentWeights", ...)}.
It is created by the \code{\link{evalWeights}} method.
}
\section{Slots}{
  \describe{
    \item{\code{w}:}{Object of class \code{"ANY"} ~~ }
    \item{\code{type}:}{Object of class \code{"character"} ~~ }
    \item{\code{wSpec}:}{Object of class \code{"list"} ~~ }
    \item{\code{Sigma}:}{Object of class \code{"ANY"} ~~ }
    \item{\code{momNames}:}{Object of class \code{"list"} ~~ }
    \item{\code{eqnNames}:}{Object of class \code{"character"} ~~ }
    \item{\code{sameMom}:}{Object of class \code{"logical"} ~~ }
  }
}
\section{Methods}{
  \describe{
    \item{[}{\code{signature(x = "sysMomentWeights", i = "missing", j = "list")}: ... }
    \item{[}{\code{signature(x = "sysMomentWeights", i = "numeric", j = "list")}: ... }
    \item{[}{\code{signature(x = "sysMomentWeights", i = "numeric", j = "missing")}: ... }
    \item{evalGmmObj}{\code{signature(object = "sysModel", theta = "list", wObj = "sysMomentWeights")}: ... }
    \item{print}{\code{signature(x = "sysMomentWeights")}: ... }
    \item{quadra}{\code{signature(w = "sysMomentWeights", x = "matrixORnumeric", y = "matrixORnumeric")}: ... }
    \item{quadra}{\code{signature(w = "sysMomentWeights", x = "matrixORnumeric", y = "missing")}: ... }
    \item{quadra}{\code{signature(w = "sysMomentWeights", x = "missing", y = "missing")}: ... }
    \item{show}{\code{signature(object = "sysMomentWeights")}: ... }
    \item{solveGmm}{\code{signature(object = "rslinearModel", wObj = "sysMomentWeights")}: ... }
    \item{solveGmm}{\code{signature(object = "slinearModel", wObj = "sysMomentWeights")}: ... }
    \item{solveGmm}{\code{signature(object = "snonlinearModel", wObj = "sysMomentWeights")}: ... }
	 }
}

\examples{
showClass("sysMomentWeights")
}
\keyword{classes}
