\name{summaryKclass-class}
\docType{class}
\alias{summaryKclass-class}
\alias{print,summaryKclass-method}
\alias{show,summaryKclass-method}

\title{Class \code{"summaryKclass"}}
\description{
The class that stores the summary statistics of model fitted by K-Class estimators.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("summaryKclass", ...)}.
It is the summary statistics of models estimated by \code{\link{kclassfit}}.
}
\section{Slots}{
  \describe{
    \item{\code{kappa}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{method}:}{Object of class \code{"character"} ~~ }
    \item{\code{origModel}:}{Object of class \code{"linearModel"} ~~ }
    \item{\code{coef}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{specTest}:}{Object of class \code{"specTest"} ~~ }
    \item{\code{strength}:}{Object of class \code{"list"} ~~ }
    \item{\code{model}:}{Object of class \code{"momentModel"} ~~ }
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
\section{Extends}{
Class \code{"\linkS4class{summaryGmm}"}, directly.
}
\section{Methods}{
  \describe{
    \item{print}{\code{signature(x = "summaryKclass")}: ... }
    \item{show}{\code{signature(object = "summaryKclass")}: ... }
	 }
}
\examples{
showClass("summaryKclass")
}
\keyword{classes}
