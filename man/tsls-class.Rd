\name{tsls-class}
\docType{class}
\alias{tsls-class}

\title{Class \code{"tsls"}}
\description{
Class that contains a fitted model using two-stage least squares
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("tsls", ...)}.
It is created my the 
}
\section{Slots}{
  \describe{
    \item{\code{theta}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{convergence}:}{Object of class \code{"numericORNULL"} ~~ }
    \item{\code{convIter}:}{Object of class \code{"numericORNULL"} ~~ }
    \item{\code{call}:}{Object of class \code{"callORNULL"} ~~ }
    \item{\code{type}:}{Object of class \code{"character"} ~~ }
    \item{\code{wObj}:}{Object of class \code{"momentWeights"} ~~ }
    \item{\code{niter}:}{Object of class \code{"integer"} ~~ }
    \item{\code{efficientGmm}:}{Object of class \code{"logical"} ~~ }
    \item{\code{model}:}{Object of class \code{"momentModel"} ~~ }
  }
}
\section{Extends}{
Class \code{"\linkS4class{gmmfit}"}, directly.
}

\examples{
showClass("tsls")
}
\keyword{classes}
