\name{kclassfit-class}
\docType{class}
\alias{kclassfit-class}
\alias{print,kclassfit-method}
\alias{show,kclassfit-method}
\alias{specTest,kclassfit,missing-method}
\alias{summary,kclassfit-method}

\title{Class \code{"kclassfit"}}
\description{
This is the object that stores the estimation result from the K-Class
estimation method. The class includes the Limited Information Maximum
Likelihood (LIML) and its modified version proposed by Fuller (1977).
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("kclassfit", ...)}.
It is created by the \code{\link{kclassfit}} function.
}
\section{Slots}{
  \describe{
    \item{\code{kappa}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{method}:}{Object of class \code{"character"} ~~ }
    \item{\code{origModel}:}{Object of class \code{"linearModel"} ~~ }
    \item{\code{theta}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{convergence}:}{Object of class \code{"numericORNULL"} ~~ }
    \item{\code{convIter}:}{Object of class \code{"numericORNULL"} ~~ }
    \item{\code{call}:}{Object of class \code{"callORNULL"} ~~ }
    \item{\code{type}:}{Object of class \code{"character"} ~~ }
    \item{\code{wObj}:}{Object of class \code{"momentWeights"} ~~ }
    \item{\code{niter}:}{Object of class \code{"integer"} ~~ }
    \item{\code{efficientGmm}:}{Object of class \code{"logical"} ~~ }
    \item{\code{model}:}{Object of class \code{"momentModel"} ~~ }
  }
}
\section{Extends}{
Class \code{"\linkS4class{gmmfit}"}, directly.
}
\section{Methods}{
  \describe{
    \item{print}{\code{signature(x = "kclassfit")}: ... }
    \item{show}{\code{signature(object = "kclassfit")}: ... }
    \item{specTest}{\code{signature(object = "kclassfit", which = "missing")}: ... }
    \item{summary}{\code{signature(object = "kclassfit")}: ... }
	 }
}


\examples{
showClass("kclassfit")
}
\keyword{classes}
