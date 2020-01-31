\name{formulaModel-class}
\docType{class}
\alias{formulaModel-class}

\title{Class \code{"formulaModel"}}
\description{
Class for moment-based models for which moments are expressed using 
formulas.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("formulaModel", ...)}.
It is generated my \code{\link{momentModel}}.
}
\section{Slots}{
  \describe{
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
Class \code{"\linkS4class{allNLModel}"}, directly.
Class \code{"\linkS4class{momentModel}"}, directly.
}
\section{Methods}{
  \describe{
    \item{[}{\code{signature(x = "formulaModel", i = "numeric", j = "missing")}: ... }
    \item{evalDMoment}{\code{signature(object = "formulaModel")}: ... }
    \item{evalMoment}{\code{signature(object = "formulaModel")}: ... }
    \item{gmmFit}{\code{signature(model = "formulaModel")}: ... }
    \item{modelDims}{\code{signature(object = "formulaModel")}: ... }
    \item{momentStrength}{\code{signature(object = "formulaModel")}: ... }
    \item{restModel}{\code{signature(object = "formulaModel")}: ... }
    \item{subset}{\code{signature(x = "formulaModel")}: ... }
	 }
}

\examples{
showClass("formulaModel")
}
\keyword{classes}
