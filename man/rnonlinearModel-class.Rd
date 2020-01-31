\name{rnonlinearModel-class}
\docType{class}
\alias{rnonlinearModel-class}

\title{Class \code{"rnonlinearModel"}}
\description{
A class for restricted moment-based models for which moment conditions
are orthogonality conditions between instruments and the residuals from
a nonlinear regression.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("rnonlinearModel",
  ...)}.
It is created by \code{\link{restModel-methods}}.
}
\section{Slots}{
  \describe{
    \item{\code{R}:}{Object of class \code{"list"} ~~ }
    \item{\code{cstSpec}:}{Object of class \code{"list"} ~~ }
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
Class \code{"\linkS4class{nonlinearModel}"}, directly.
Class \code{"\linkS4class{rmomentModel}"}, directly.
Class \code{"\linkS4class{regModel}"}, by class "nonlinearModel", distance 2.
Class \code{"\linkS4class{allNLModel}"}, by class "nonlinearModel", distance 2.
Class \code{"\linkS4class{momentModel}"}, by class "nonlinearModel", distance 2.
}
\section{Methods}{
  \describe{
    \item{coef}{\code{signature(object = "rnonlinearModel")}: ... }
    \item{evalDMoment}{\code{signature(object = "rnonlinearModel")}: ... }
    \item{getRestrict}{\code{signature(object = "rnonlinearModel")}: ... }
    \item{gmmFit}{\code{signature(model = "rnonlinearModel")}: ... }
    \item{modelDims}{\code{signature(object = "rnonlinearModel")}: ... }
    \item{print}{\code{signature(x = "rnonlinearModel")}: ... }
    \item{printRestrict}{\code{signature(object = "rnonlinearModel")}: ... }
	 }
}

\examples{
showClass("rnonlinearModel")
}
\keyword{classes}
