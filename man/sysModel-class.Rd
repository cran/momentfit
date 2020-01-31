\name{sysModel-class}
\docType{class}
\alias{sysModel-class}

\title{Class \code{"sysModel"}}
\description{
A union class for all systems of equations.
}
\section{Objects from the Class}{A virtual Class: No objects may be created from it.}
\section{Methods}{
  \describe{
    \item{[}{\code{signature(x = "sysModel", i = "missing", j = "list")}: ... }
    \item{[}{\code{signature(x = "sysModel", i = "missing", j = "missing")}: ... }
    \item{[}{\code{signature(x = "sysModel", i = "numeric", j = "list")}: ... }
    \item{Dresiduals}{\code{signature(object = "sysModel")}: ... }
    \item{evalDMoment}{\code{signature(object = "sysModel")}: ... }
    \item{evalGmmObj}{\code{signature(object = "sysModel", theta = "list", wObj = "sysMomentWeights")}: ... }
    \item{evalMoment}{\code{signature(object = "sysModel")}: ... }
    \item{evalWeights}{\code{signature(object = "sysModel")}: ... }
    \item{getRestrict}{\code{signature(object = "sysModel")}: ... }
    \item{gmmFit}{\code{signature(model = "sysModel")}: ... }
    \item{print}{\code{signature(x = "sysModel")}: ... }
    \item{residuals}{\code{signature(object = "sysModel")}: ... }
    \item{show}{\code{signature(object = "sysModel")}: ... }
    \item{subset}{\code{signature(x = "sysModel")}: ... }
    \item{vcov}{\code{signature(object = "sysModel")}: ... }
	 }
}
\examples{
showClass("sysModel")
}
\keyword{classes}
