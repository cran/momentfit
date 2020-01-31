\name{rsnonlinearModel-class}
\docType{class}
\alias{rsnonlinearModel-class}

\title{Class \code{"rsnonlinearModel"}}
\description{
A class for restricted systems of nonlinear equations.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("rsnonlinearModel", ...)}.
It is created by \code{\link{restModel-methods}}.
}
\section{Slots}{
  \describe{
    \item{\code{R}:}{Object of class \code{"list"} ~~ }
    \item{\code{cstSpec}:}{Object of class \code{"list"} ~~ }
    \item{\code{data}:}{Object of class \code{"data.frame"} ~~ }
    \item{\code{instT}:}{Object of class \code{"list"} ~~ }
    \item{\code{vcov}:}{Object of class \code{"character"} ~~ }
    \item{\code{theta0}:}{Object of class \code{"list"} ~~ }
    \item{\code{n}:}{Object of class \code{"integer"} ~~ }
    \item{\code{q}:}{Object of class \code{"integer"} ~~ }
    \item{\code{k}:}{Object of class \code{"integer"} ~~ }
    \item{\code{parNames}:}{Object of class \code{"list"} ~~ }
    \item{\code{momNames}:}{Object of class \code{"list"} ~~ }
    \item{\code{fRHS}:}{Object of class \code{"list"} ~~ }
    \item{\code{fLHS}:}{Object of class \code{"list"} ~~ }
    \item{\code{eqnNames}:}{Object of class \code{"character"} ~~ }
    \item{\code{vcovOptions}:}{Object of class \code{"list"} ~~ }
    \item{\code{centeredVcov}:}{Object of class \code{"logical"} ~~ }
    \item{\code{sameMom}:}{Object of class \code{"logical"} ~~ }
    \item{\code{SUR}:}{Object of class \code{"logical"} ~~ }
    \item{\code{varNames}:}{Object of class \code{"list"} ~~ }
    \item{\code{isEndo}:}{Object of class \code{"list"} ~~ }
    \item{\code{omit}:}{Object of class \code{"integer"} ~~ }
    \item{\code{survOptions}:}{Object of class \code{"list"} ~~ }
    \item{\code{sSpec}:}{Object of class \code{"sSpec"} ~~ }
    \item{\code{smooth}:}{Object of class \code{"logical"} ~~ }
  }
}
\section{Extends}{
Class \code{"\linkS4class{snonlinearModel}"}, directly.
Class \code{"\linkS4class{rsysModel}"}, directly.
Class \code{"\linkS4class{sysModel}"}, by class "snonlinearModel", distance 2.
}
\section{Methods}{
No methods defined with class "rsnonlinearModel" in the signature.
}

\examples{
showClass("rsnonlinearModel")
}
\keyword{classes}
