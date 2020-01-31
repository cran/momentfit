\name{hypothesisTest-class}
\docType{class}
\alias{hypothesisTest-class}

\title{Class \code{"hypothesisTest"}}
\description{
A class to store results form an hypothesis test. 
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("hypothesisTest", ...)}.
It is created by \code{\link{hypothesisTest}}.
}
\section{Slots}{
  \describe{
    \item{\code{test}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{hypothesis}:}{Object of class \code{"character"} ~~ }
    \item{\code{dist}:}{Object of class \code{"character"} ~~ }
    \item{\code{df}:}{Object of class \code{"integer"} ~~ }
    \item{\code{pvalue}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{type}:}{Object of class \code{"character"} ~~ }
  }
}
\section{Methods}{
  \describe{
    \item{print}{\code{signature(x = "hypothesisTest")}: ... }
    \item{show}{\code{signature(object = "hypothesisTest")}: ... }
	 }
}

\examples{
showClass("hypothesisTest")
}
\keyword{classes}
