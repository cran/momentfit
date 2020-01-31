\name{nonlinearModel-class}
\Rdversion{1.1}
\docType{class}
\alias{nonlinearModel-class}

\title{Class \code{"nonlinearModel"}}
\description{
  Class for moment-based models for which moment conditions are
  orthogonality conditions between instruments and the residuals from a
  nonlinear regression. 
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("nonlinearModel",
  ...)}.
It is generated my \code{\link{momentModel}}.
}
\section{Slots}{
  \describe{
    \item{\code{modelF}:}{Object of class \code{"data.frame"} ~~ }
    \item{\code{instF}:}{Object of class \code{"data.frame"} ~~ }
    \item{\code{vcov}:}{Object of class \code{"character"} ~~ }
    \item{\code{theta0}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{n}:}{Object of class \code{"integer"} ~~ }
    \item{\code{q}:}{Object of class \code{"integer"} ~~ }
    \item{\code{k}:}{Object of class \code{"integer"} ~~ }
    \item{\code{parNames}:}{Object of class \code{"character"} ~~ }
    \item{\code{momNames}:}{Object of class \code{"character"} ~~ }
    \item{\code{fRHS}:}{Object of class \code{"expression"} ~~ }
    \item{\code{fLHS}:}{Object of class \code{"expressionORNULL"} ~~ }
    \item{\code{vcovOptions}:}{Object of class \code{"list"} ~~ }
    \item{\code{centeredVcov}:}{Object of class \code{"logical"} ~~ }
    \item{\code{varNames}:}{Object of class \code{"character"} ~~ }
    \item{\code{isEndo}:}{Object of class \code{"logical"} ~~ }
    \item{\code{omit}:}{Object of class \code{"integer"} ~~ }
    \item{\code{survOptions}:}{Object of class \code{"list"} ~~ }
    \item{\code{sSpec}:}{Object of class \code{"sSpec"} ~~ }
    \item{\code{smooth}:}{Object of class \code{"logical"} ~~ }
  }
}
\section{Extends}{
Class \code{"\linkS4class{regModel}"}, directly.
Class \code{"\linkS4class{allNLModel}"}, directly.
Class \code{"\linkS4class{momentModel}"}, directly.
}
\section{Methods}{
  \describe{
    \item{Dresiduals}{\code{signature(object = "nonlinearModel")}: ... }
    \item{merge}{\code{signature(x = "nonlinearModel", y = "nonlinearModel")}: ... }
    \item{merge}{\code{signature(x = "snonlinearModel", y = "nonlinearModel")}: ... }
    \item{model.matrix}{\code{signature(object = "nonlinearModel")}: ... }
    \item{modelDims}{\code{signature(object = "nonlinearModel")}: ... }
    \item{momentStrength}{\code{signature(object = "nonlinearModel")}: ... }
    \item{residuals}{\code{signature(object = "nonlinearModel")}: ... }
    \item{restModel}{\code{signature(object = "nonlinearModel")}: ... }
	 }
}

\examples{
showClass("nonlinearModel")
}
\keyword{classes}
