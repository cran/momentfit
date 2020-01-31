\name{regModel-class}
\docType{class}
\alias{regModel-class}

\title{Class \code{"regModel"}}
\description{
A union class for \code{"linearModel"} and \code{"nonlinearModel"}
classes.
}
\section{Objects from the Class}{A virtual Class: No objects may be created from it.}
\section{Methods}{
  \describe{
    \item{[}{\code{signature(x = "regModel", i = "numeric", j = "missing")}: ... }
    \item{evalDMoment}{\code{signature(object = "regModel")}: ... }
    \item{evalMoment}{\code{signature(object = "regModel")}: ... }
    \item{subset}{\code{signature(x = "regModel")}: ... }
	 }
}

\examples{
showClass("regModel")
}
\keyword{classes}
