\name{rfunctionModel-class}
\docType{class}
\alias{rfunctionModel-class}

\title{Class \code{"rfunctionModel"}}
\description{
A restricted moment-based model for which moment conditions are defined
by a function. 
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("rfunctionModel", ...)}.
It is created by \code{\link{restModel-methods}}. 
}
\section{Slots}{
  \describe{
    \item{\code{R}:}{Object of class \code{"list"} ~~ }
    \item{\code{cstSpec}:}{Object of class \code{"list"} ~~ }
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
Class \code{"\linkS4class{functionModel}"}, directly.
Class \code{"\linkS4class{rmomentModel}"}, directly.
Class \code{"\linkS4class{allNLModel}"}, by class "functionModel", distance 2.
Class \code{"\linkS4class{momentModel}"}, by class "functionModel", distance 2.
}
\section{Methods}{
  \describe{
    \item{[}{\code{signature(x = "rfunctionModel", i = "numeric", j = "missing")}: ... }
    \item{coef}{\code{signature(object = "rfunctionModel")}: ... }
    \item{evalDMoment}{\code{signature(object = "rfunctionModel")}: ... }
    \item{getRestrict}{\code{signature(object = "rfunctionModel")}: ... }
    \item{modelDims}{\code{signature(object = "rfunctionModel")}: ... }
    \item{print}{\code{signature(x = "rfunctionModel")}: ... }
    \item{printRestrict}{\code{signature(object = "rfunctionModel")}: ... }
	 }
}

\examples{
showClass("rfunctionModel")
}
\keyword{classes}
