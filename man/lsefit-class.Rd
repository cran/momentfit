\name{lsefit-class}
\Rdversion{1.1}
\docType{class}
\alias{lsefit-class}
\alias{print,lsefit-method}
\alias{show,lsefit-method}

\title{Class \code{"lsefit"}}
\description{
A class for least squares estimate of different \code{momentModel}
objects. 
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("lsefit", ...)}.
It is created by \code{\link{lse}}. It includes the information about
the model being estimated and the estimation based on \code{\link{lm}}.
}
\section{Slots}{
  \describe{
    \item{\code{model}:}{Object of class \code{"linearModel"} ~~ }
    \item{\code{.S3Class}:}{Object of class \code{"character"} ~~ }
  }
}
\section{Extends}{
Class \code{"\linkS4class{lm}"}, directly.
Class \code{"\linkS4class{oldClass}"}, by class "lm", distance 2.
}
\section{Methods}{
  \describe{
    \item{print}{\code{signature(x = "lsefit")}: ... }
    \item{show}{\code{signature(object = "lsefit")}: ... }
	 }
}

\examples{
showClass("lsefit")
}
\keyword{classes}
