\name{functionModel-class}
\docType{class}
\alias{functionModel-class}

\title{Class \code{"functionModel"}}
\description{
Class for moment-based models for which moment conditions are defined
using a function.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("functionModel", ...)}.
It is generated my \code{\link{momentModel}}.
}
\section{Slots}{
  \describe{
    \item{\code{X}:}{Object of class \code{"ANY"} ~~ }
    \item{\code{fct}:}{Object of class \code{"function"} ~~ }
    \item{\code{dfct}:}{Object of class \code{"functionORNULL"} ~~ }
    \item{\code{vcov}:}{Object of class \code{"character"} ~~ }
    \item{\code{theta0}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{n}:}{Object of class \code{"integer"} ~~ }
    \item{\code{q}:}{Object of class \code{"integer"} ~~ }
    \item{\code{k}:}{Object of class \code{"integer"} ~~ }
    \item{\code{parNames}:}{Object of class \code{"character"} ~~ }
    \item{\code{momNames}:}{Object of class \code{"character"} ~~ }
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
Class \code{"\linkS4class{allNLModel}"}, directly.
Class \code{"\linkS4class{momentModel}"}, directly.
}
\section{Methods}{
  \describe{
    \item{[}{\code{signature(x = "functionModel", i = "numeric", j = "missing")}: ... }
    \item{evalDMoment}{\code{signature(object = "functionModel")}: ... }
    \item{evalMoment}{\code{signature(object = "functionModel")}: ... }
    \item{modelDims}{\code{signature(object = "functionModel")}: ... }
    \item{momentStrength}{\code{signature(object = "functionModel")}: ... }
    \item{restModel}{\code{signature(object = "functionModel")}: ... }
    \item{subset}{\code{signature(x = "functionModel")}: ... }
	 }
}

\examples{
showClass("functionModel")
}
\keyword{classes}
