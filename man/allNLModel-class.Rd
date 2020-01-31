\name{allNLModel-class}
\docType{class}
\alias{allNLModel-class}


\title{Class \code{"allNLModel"}}
\description{
  A union class for all nonlinear models. It includes
  \code{"nonlinearModel"}, \code{"formulaModel"}, and
  \code{"functionModel"}. 
}
\section{Objects from the Class}{A virtual Class: No objects may be created from it.}
\section{Methods}{
  \describe{
    \item{solveGmm}{\code{signature(object = "allNLModel", wObj = "momentWeights")}: ... }
	 }
}

\examples{
showClass("allNLModel")
}
\keyword{classes}
