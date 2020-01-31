\name{stsls-class}
\docType{class}
\alias{stsls-class}

\title{Class \code{"stsls"}}
\description{
A class to store a fitted system of equations obtained using the
two-stage least squares method.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("stsls", ...)}.
It is created my \code{\link{tsls-methods}}.
}
\section{Slots}{
  \describe{
    \item{\code{theta}:}{Object of class \code{"list"} ~~ }
    \item{\code{convergence}:}{Object of class \code{"numericORNULL"} ~~ }
    \item{\code{convIter}:}{Object of class \code{"numericORNULL"} ~~ }
    \item{\code{call}:}{Object of class \code{"callORNULL"} ~~ }
    \item{\code{type}:}{Object of class \code{"character"} ~~ }
    \item{\code{wObj}:}{Object of class \code{"sysMomentWeights"} ~~ }
    \item{\code{niter}:}{Object of class \code{"integer"} ~~ }
    \item{\code{efficientGmm}:}{Object of class \code{"logical"} ~~ }
    \item{\code{model}:}{Object of class \code{"sysModel"} ~~ }
  }
}
\section{Extends}{
Class \code{"\linkS4class{sgmmfit}"}, directly.
}
\section{Methods}{
No methods defined with class "stsls" in the signature.
}

\examples{
showClass("stsls")
}
\keyword{classes}
