\name{rformulaModel-class}
\docType{class}
\alias{rformulaModel-class}

\title{Class \code{"rformulaModel"}}
\description{
A class for restricted moment-based models for which moment conditions
are expressed using a list of formulas.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("rformulaModel", ...)}.
It is created by \code{\link{restModel-methods}}.
}
\section{Slots}{
  \describe{
    \item{\code{R}:}{Object of class \code{"list"} ~~ }
    \item{\code{cstSpec}:}{Object of class \code{"list"} ~~ }
    \item{\code{modelF}:}{Object of class \code{"data.frame"} ~~ }
    \item{\code{vcov}:}{Object of class \code{"character"} ~~ }
    \item{\code{theta0}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{n}:}{Object of class \code{"integer"} ~~ }
    \item{\code{q}:}{Object of class \code{"integer"} ~~ }
    \item{\code{k}:}{Object of class \code{"integer"} ~~ }
    \item{\code{parNames}:}{Object of class \code{"character"} ~~ }
    \item{\code{momNames}:}{Object of class \code{"character"} ~~ }
    \item{\code{fRHS}:}{Object of class \code{"list"} ~~ }
    \item{\code{fLHS}:}{Object of class \code{"list"} ~~ }
    \item{\code{vcovOptions}:}{Object of class \code{"list"} ~~ }
    \item{\code{centeredVcov}:}{Object of class \code{"logical"} ~~ }
    \item{\code{varNames}:}{Object of class \code{"character"} ~~ }
    \item{\code{isEndo}:}{Object of class \code{"logical"} ~~ }
    \item{\code{isMDE}:}{Object of class \code{"logical"} ~~ }
    \item{\code{omit}:}{Object of class \code{"integer"} ~~ }
    \item{\code{survOptions}:}{Object of class \code{"list"} ~~ }
    \item{\code{sSpec}:}{Object of class \code{"sSpec"} ~~ }
    \item{\code{smooth}:}{Object of class \code{"logical"} ~~ }
  }
}
\section{Extends}{
Class \code{"\linkS4class{formulaModel}"}, directly.
Class \code{"\linkS4class{rmomentModel}"}, directly.
Class \code{"\linkS4class{allNLModel}"}, by class "formulaModel", distance 2.
Class \code{"\linkS4class{momentModel}"}, by class "formulaModel", distance 2.
}
\section{Methods}{
  \describe{
    \item{coef}{\code{signature(object = "rformulaModel")}: ... }
    \item{evalDMoment}{\code{signature(object = "rformulaModel")}: ... }
    \item{getRestrict}{\code{signature(object = "rformulaModel")}: ... }
    \item{gmmFit}{\code{signature(model = "rformulaModel")}: ... }
    \item{modelDims}{\code{signature(object = "rformulaModel")}: ... }
    \item{print}{\code{signature(x = "rformulaModel")}: ... }
    \item{printRestrict}{\code{signature(object = "rformulaModel")}: ... }
	 }
}

\examples{
showClass("rformulaModel")
}
\keyword{classes}
