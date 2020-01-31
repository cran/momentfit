\name{mconfint-class}
\docType{class}
\alias{mconfint-class}

\title{Class \code{"mconfint"}}
\description{
A class to store confidence region.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("mconfint", ...)}.
It is created by the \code{"confint"} method with the option
\code{area=TRUE} (see \code{\link{confint-methods}}).
}
\section{Slots}{
  \describe{
    \item{\code{areaPoints}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{type}:}{Object of class \code{"character"} ~~ }
    \item{\code{level}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{theta}:}{Object of class \code{"numeric"} ~~ }
  }
}
\section{Methods}{
  \describe{
    \item{plot}{\code{signature(x = "mconfint")}: ... }
    \item{print}{\code{signature(x = "mconfint")}: ... }
    \item{show}{\code{signature(object = "mconfint")}: ... }
	 }
}

\examples{
showClass("mconfint")
}
\keyword{classes}
