\name{specTest-class}
\docType{class}
\alias{specTest-class}

\title{Class \code{"specTest"}}
\description{
A class to store results from a specification test.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("specTest", ...)}.
It is created my \code{\link{specTest-methods}}.
}
\section{Slots}{
  \describe{
    \item{\code{test}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{testname}:}{Object of class \code{"character"} ~~ }
  }
}
\section{Methods}{
  \describe{
    \item{print}{\code{signature(x = "specTest")}: ... }
    \item{show}{\code{signature(object = "specTest")}: ... }
	 }
}

\examples{
showClass("specTest")
}
\keyword{classes}
