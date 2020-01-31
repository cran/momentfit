\name{linearModel-class}
\docType{class}
\alias{linearModel-class}

\title{Class \code{"linearModel"}}
\description{
Class for moment-based models for which moment conditions are linear and
expressed by a formula.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("linearModel",
  ...)}.
It is generated my \code{\link{momentModel}}.
}
\section{Slots}{
  \describe{
    \item{\code{modelF}:}{Object of class \code{"data.frame"} ~~ }
    \item{\code{instF}:}{Object of class \code{"data.frame"} ~~ }
    \item{\code{vcov}:}{Object of class \code{"character"} ~~ }
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
Class \code{"\linkS4class{regModel}"}, directly.
Class \code{"\linkS4class{momentModel}"}, directly.
}
\section{Methods}{
  \describe{
    \item{Dresiduals}{\code{signature(object = "linearModel")}: ... }
    \item{merge}{\code{signature(x = "linearModel", y = "linearModel")}: ... }
    \item{merge}{\code{signature(x = "slinearModel", y = "linearModel")}: ... }
    \item{model.matrix}{\code{signature(object = "linearModel")}: ... }
    \item{modelDims}{\code{signature(object = "linearModel")}: ... }
    \item{modelResponse}{\code{signature(object = "linearModel")}: ... }
    \item{momentStrength}{\code{signature(object = "linearModel")}: ... }
    \item{residuals}{\code{signature(object = "linearModel")}: ... }
    \item{restModel}{\code{signature(object = "linearModel")}: ... }
    \item{solveGmm}{\code{signature(object = "linearModel", wObj = "momentWeights")}: ... }
    \item{tsls}{\code{signature(model = "linearModel")}: ... }
	 }
}

\examples{
showClass("linearModel")
}
\keyword{classes}
