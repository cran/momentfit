\name{rslinearModel-class}
\docType{class}
\alias{rslinearModel-class}


\title{Class \code{"rslinearModel"}}
\description{
A class for restricted system of linear equations.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("rslinearModel",
  ...)}.
It is created by \code{\link{restModel-methods}}.
}
\section{Slots}{
  \describe{
    \item{\code{cstLHS}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{cstRHS}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{cstSpec}:}{Object of class \code{"list"} ~~ }
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
Class \code{"\linkS4class{slinearModel}"}, directly.
Class \code{"\linkS4class{rsysModel}"}, directly.
Class \code{"\linkS4class{sysModel}"}, by class "slinearModel", distance 2.
}
\section{Methods}{
  \describe{
    \item{[}{\code{signature(x = "rslinearModel", i = "numeric", j = "missing")}: ... }
    \item{coef}{\code{signature(object = "rslinearModel")}: ... }
    \item{evalDMoment}{\code{signature(object = "rslinearModel")}: ... }
    \item{evalMoment}{\code{signature(object = "rslinearModel")}: ... }
    \item{evalWeights}{\code{signature(object = "rslinearModel")}: ... }
    \item{getRestrict}{\code{signature(object = "rslinearModel")}: ... }
    \item{gmmFit}{\code{signature(model = "rslinearModel")}: ... }
    \item{model.matrix}{\code{signature(object = "rslinearModel")}: ... }
    \item{modelDims}{\code{signature(object = "rslinearModel")}: ... }
    \item{modelResponse}{\code{signature(object = "rslinearModel")}: ... }
    \item{print}{\code{signature(x = "rslinearModel")}: ... }
    \item{printRestrict}{\code{signature(object = "rslinearModel")}: ... }
    \item{residuals}{\code{signature(object = "rslinearModel")}: ... }
    \item{solveGmm}{\code{signature(object = "rslinearModel", wObj = "sysMomentWeights")}: ... }
    \item{ThreeSLS}{\code{signature(model = "rslinearModel")}: ... }
	 }
}

\examples{
showClass("rslinearModel")
}
\keyword{classes}
