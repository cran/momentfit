\name{gelfit-class}
\docType{class}
\alias{gelfit-class}

\title{Class \code{"gelfit"}}
\description{
A class to store fitted models obtained using a GEL method.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("gelfit", ...)}.
It is created by \code{\link{gelFit}}.
}
\section{Slots}{
  \describe{
    \item{\code{theta}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{convergence}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{lambda}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{lconvergence}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{call}:}{Object of class \code{"callORNULL"} ~~ }
    \item{\code{gelType}:}{Object of class \code{"list"} ~~ }
    \item{\code{vcov}:}{Object of class \code{"list"} ~~ }
    \item{\code{model}:}{Object of class \code{"momentModel"} ~~ }
    \item{\code{restrictedLam}:}{Object of class \code{"integer"} ~~ }
  }
}
\section{Methods}{
  \describe{
    \item{coef}{\code{signature(object = "gelfit")}: ... }
    \item{confint}{\code{signature(object = "gelfit")}: ... }
    \item{getImpProb}{\code{signature(object = "gelfit")}: ... }
    \item{momFct}{\code{signature(eta = "numeric", object = "gelfit")}: ... }
    \item{print}{\code{signature(x = "gelfit")}: ... }
    \item{residuals}{\code{signature(object = "gelfit")}: ... }
    \item{show}{\code{signature(object = "gelfit")}: ... }
    \item{specTest}{\code{signature(object = "gelfit", which = "missing")}: ... }
    \item{summary}{\code{signature(object = "gelfit")}: ... }
    \item{update}{\code{signature(object = "gelfit")}: ... }
    \item{vcov}{\code{signature(object = "gelfit")}: ... }
	 }
}

\examples{
showClass("gelfit")
}
\keyword{classes}
