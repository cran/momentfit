\name{summarySysGmm-class}
\docType{class}
\alias{summarySysGmm-class}

\title{Class \code{"summarySysGmm"}}
\description{
A class to store the summary of a system of equations fitted by GMM.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("summarySysGmm", ...)}.
It is created by \code{\link{summary-methods}}.
}
\section{Slots}{
  \describe{
    \item{\code{coef}:}{Object of class \code{"list"} ~~ }
    \item{\code{specTest}:}{Object of class \code{"specTest"} ~~ }
    \item{\code{strength}:}{Object of class \code{"list"} ~~ }
    \item{\code{model}:}{Object of class \code{"sysModel"} ~~ }
    \item{\code{sandwich}:}{Object of class \code{"logical"} ~~ }
    \item{\code{type}:}{Object of class \code{"character"} ~~ }
    \item{\code{convergence}:}{Object of class \code{"numericORNULL"} ~~ }
    \item{\code{convIter}:}{Object of class \code{"numericORNULL"} ~~ }
    \item{\code{wSpec}:}{Object of class \code{"list"} ~~ }
    \item{\code{niter}:}{Object of class \code{"integer"} ~~ }
    \item{\code{df.adj}:}{Object of class \code{"logical"} ~~ }
    \item{\code{breadOnly}:}{Object of class \code{"logical"} ~~ }
  }
}
\section{Methods}{
  \describe{
    \item{print}{\code{signature(x = "summarySysGmm")}: ... }
    \item{show}{\code{signature(object = "summarySysGmm")}: ... }
	 }
}

\examples{
showClass("summarySysGmm")
}
\keyword{classes}
