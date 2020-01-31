\name{sgmmfit-class}
\docType{class}
\alias{sgmmfit-class}

\title{Class \code{"sgmmfit"}}
\description{
  Class to store fitted system of equations obtained using the GMM method.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("sgmmfit", ...)}.
It is created by \code{\link{gmmFit}}.
}
\section{Slots}{
  \describe{
    \item{\code{theta}:}{Object of class \code{"list"} ~~ }
    \item{\code{convergence}:}{Object of class \code{"numericORNULL"} ~~ }
    \item{\code{convIter}:}{Object of class \code{"numericORNULL"} ~~ }
    \item{\code{call}:}{Object of class \code{"callORNULL"} ~~ }
    \item{\code{type}:}{Object of class \code{"character"} ~~ }
    \item{\code{wObj}:}{Object of class \code{"sysMomentWeights"} ~~ }
    \item{\code{niter}:}{Object of class \code{"integer"} ~~ }
    \item{\code{efficientGmm}:}{Object of class \code{"logical"} ~~ }
    \item{\code{model}:}{Object of class \code{"sysModel"} ~~ }
  }
}
\section{Methods}{
  \describe{
    \item{bread}{\code{signature(x = "sgmmfit")}: ... }
    \item{coef}{\code{signature(object = "sgmmfit")}: ... }
    \item{hypothesisTest}{\code{signature(object.u = "missing", object.r = "sgmmfit")}: ... }
    \item{hypothesisTest}{\code{signature(object.u = "sgmmfit", object.r = "missing")}: ... }
    \item{hypothesisTest}{\code{signature(object.u = "sgmmfit", object.r = "sgmmfit")}: ... }
    \item{meatGmm}{\code{signature(object = "sgmmfit")}: ... }
    \item{print}{\code{signature(x = "sgmmfit")}: ... }
    \item{residuals}{\code{signature(object = "sgmmfit")}: ... }
    \item{show}{\code{signature(object = "sgmmfit")}: ... }
    \item{specTest}{\code{signature(object = "sgmmfit", which = "missing")}: ... }
    \item{summary}{\code{signature(object = "sgmmfit")}: ... }
    \item{vcov}{\code{signature(object = "sgmmfit")}: ... }
	 }
}

\examples{
showClass("sgmmfit")
}
\keyword{classes}
