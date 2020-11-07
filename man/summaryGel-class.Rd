\name{summaryGel-class}
\docType{class}
\alias{summaryGel-class}

\title{Class \code{"summaryGel"}}
\description{
Class to store the summary of a model fitted by GEL.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("summaryGel", ...)}.
It is created by \code{link{summary-methods}}.
}
\section{Slots}{
  \describe{
    \item{\code{coef}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{specTest}:}{Object of class \code{"specTest"} ~~ }
    \item{\code{model}:}{Object of class \code{"momentModel"} ~~ }
    \item{\code{lambda}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{convergence}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{lconvergence}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{impProb}:}{Object of class \code{"list"} ~~ }
    \item{\code{gelType}:}{Object of class \code{"list"} ~~ }
    \item{\code{restrictedLam}:}{Object of class \code{"integer"} ~~ }
  }
}
\section{Methods}{
  \describe{
    \item{print}{\code{signature(x = "summaryGel")}: ... }
    \item{show}{\code{signature(object = "summaryGel")}: ... }
	 }
}

\examples{
showClass("summaryGel")
}
\keyword{classes}
