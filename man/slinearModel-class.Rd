\name{slinearModel-class}
\docType{class}
\alias{slinearModel-class}

\title{Class \code{"slinearModel"}}
\description{
A class for systems of linear equations.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("slinearModel", ...)}.
It is created by \code{\link{momentModel}}.
}
\section{Slots}{
  \describe{
    \item{\code{modelT}:}{Object of class \code{"list"} ~~ }
    \item{\code{instT}:}{Object of class \code{"list"} ~~ }
    \item{\code{data}:}{Object of class \code{"data.frame"} ~~ }
    \item{\code{vcov}:}{Object of class \code{"character"} ~~ }
    \item{\code{n}:}{Object of class \code{"integer"} ~~ }
    \item{\code{q}:}{Object of class \code{"integer"} ~~ }
    \item{\code{k}:}{Object of class \code{"integer"} ~~ }
    \item{\code{parNames}:}{Object of class \code{"list"} ~~ }
    \item{\code{momNames}:}{Object of class \code{"list"} ~~ }
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
Class \code{"\linkS4class{sysModel}"}, directly.
}
\section{Methods}{
  \describe{
    \item{[}{\code{signature(x = "slinearModel", i = "numeric", j = "missing")}: ... }
    \item{merge}{\code{signature(x = "slinearModel", y = "linearModel")}: ... }
    \item{model.matrix}{\code{signature(object = "slinearModel")}: ... }
    \item{modelDims}{\code{signature(object = "slinearModel")}: ... }
    \item{modelResponse}{\code{signature(object = "slinearModel")}: ... }
    \item{restModel}{\code{signature(object = "slinearModel")}: ... }
    \item{solveGmm}{\code{signature(object = "slinearModel", wObj = "sysMomentWeights")}: ... }
    \item{ThreeSLS}{\code{signature(model = "slinearModel")}: ... }
    \item{tsls}{\code{signature(model = "slinearModel")}: ... }
	 }
}

\examples{
showClass("slinearModel")
}
\keyword{classes}
