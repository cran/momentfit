\name{rlinearModel-class}
\Rdversion{1.1}
\docType{class}
\alias{rlinearModel-class}

\title{Class \code{"rlinearModel"}}
\description{
A class for restricted moment-based models for which moment conditions
are orthogonality conditions between instruments and the residuals from
a linear regression. 
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("rlinearModel", ...)}.  
It is created by \code{\link{restModel-methods}}.
}

\section{Slots}{
  \describe{
    \item{\code{cstLHS}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{cstRHS}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{cstSpec}:}{Object of class \code{"list"} ~~ }
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
Class \code{"\linkS4class{linearModel}"}, directly.
Class \code{"\linkS4class{rmomentModel}"}, directly.
Class \code{"\linkS4class{regModel}"}, by class "linearModel", distance 2.
Class \code{"\linkS4class{momentModel}"}, by class "linearModel", distance 2.
}
\section{Methods}{
  \describe{
    \item{coef}{\code{signature(object = "rlinearModel")}: ... }
    \item{getRestrict}{\code{signature(object = "rlinearModel")}: ... }
    \item{gmmFit}{\code{signature(model = "rlinearModel")}: ... }
    \item{model.matrix}{\code{signature(object = "rlinearModel")}: ... }
    \item{modelDims}{\code{signature(object = "rlinearModel")}: ... }
    \item{modelResponse}{\code{signature(object = "rlinearModel")}: ... }
    \item{momentStrength}{\code{signature(object = "rlinearModel")}: ... }
    \item{print}{\code{signature(x = "rlinearModel")}: ... }
    \item{printRestrict}{\code{signature(object = "rlinearModel")}: ... }
	 }
}

\examples{
showClass("rlinearModel")
}
\keyword{classes}
