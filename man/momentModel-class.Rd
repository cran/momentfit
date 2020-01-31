\name{momentModel-class}
\docType{class}
\alias{momentModel-class}

\title{Class \code{"momentModel"}}
\description{
A union class for all moment based models. It is created by
\code{\link{momentModel}}.
}
\section{Objects from the Class}{A virtual Class: No objects may be created from it.}
\section{Methods}{
  \describe{
    \item{[}{\code{signature(x = "momentModel", i = "missing", j = "missing")}: ... }
    \item{coef}{\code{signature(object = "momentModel")}: ... }
    \item{evalGel}{\code{signature(model = "momentModel")}: ... }
    \item{evalGelObj}{\code{signature(object = "momentModel", theta = "numeric", lambda = "numeric")}: ... }
    \item{evalGmm}{\code{signature(model = "momentModel")}: ... }
    \item{evalGmmObj}{\code{signature(object = "momentModel", theta = "numeric", wObj = "momentWeights")}: ... }
    \item{evalWeights}{\code{signature(object = "momentModel")}: ... }
    \item{gelFit}{\code{signature(model = "momentModel")}: ... }
    \item{getRestrict}{\code{signature(object = "momentModel")}: ... }
    \item{gmmFit}{\code{signature(model = "momentModel")}: ... }
    \item{kernapply}{\code{signature(x = "momentModel")}: ... }
    \item{print}{\code{signature(x = "momentModel")}: ... }
    \item{show}{\code{signature(object = "momentModel")}: ... }
    \item{solveGel}{\code{signature(object = "momentModel")}: ... }
    \item{update}{\code{signature(object = "momentModel")}: ... }
    \item{vcov}{\code{signature(object = "momentModel")}: ... }
    \item{vcovHAC}{\code{signature(x = "momentModel")}: ... }
	 }
}

\examples{
showClass("momentModel")
}
\keyword{classes}
