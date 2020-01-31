\name{sSpec-class}
\docType{class}
\alias{sSpec-class}

\title{Class \code{"sSpec"}}
\description{
A class to store the specifications of the kernel used to smooth moment
conditions. 
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("sSpec", ...)}.
It is created by \code{\link{kernapply-methods}}.
}
\section{Slots}{
  \describe{
    \item{\code{k}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{kernel}:}{Object of class \code{"character"} ~~ }
    \item{\code{bw}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{w}:}{Object of class \code{"tskernel"} ~~ }
    \item{\code{bwMet}:}{Object of class \code{"character"} ~~ }
  }
}
\section{Methods}{
  \describe{
    \item{print}{\code{signature(x = "sSpec")}: ... }
    \item{show}{\code{signature(object = "sSpec")}: ... }
	 }
}

\examples{
showClass("sSpec")
}
\keyword{classes}
