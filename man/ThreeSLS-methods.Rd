\name{ThreeSLS-methods}
\docType{methods}
\alias{ThreeSLS}
\alias{ThreeSLS-methods}
\alias{ThreeSLS,slinearModel-method}
\alias{ThreeSLS,rslinearModel-method}
\title{ ~~ Methods for Function \code{ThreeSLS} in Package \pkg{momentfit} ~~}
\description{
Method to estimate system of equations by Three-Stage least squares
(3SLS) or, as a special case, by Seemingly Unrelatd Regressions (SUR).
}
\usage{
\S4method{ThreeSLS}{slinearModel}(model, coefOnly=FALSE, qrZ=NULL,
Sigma=NULL)
\S4method{ThreeSLS}{rslinearModel}(model, coefOnly=FALSE, qrZ=NULL,
Sigma=NULL)
}
\arguments{
  \item{model}{An object of class \code{"slinearModel"} in which
    instruments are the same in each equation and the error terms are
    homoscedastic.}
  \item{coefOnly}{Should the method return the only the coefficients or
    create an object of class \code{"sgmmfit"}.}
  \item{qrZ}{The qr decomposition of the common instruments. It is
    mostly used by \code{\link{gmmFit}} to avoid recomputing it in
    iterative GMM or CUE. It should not be used directly unless the user
    knows what he is doing.}
  \item{Sigma}{The covariance matrix of the residuals. If not provided,
    it is computed using the residuals of the equation by equation
    two-stage least squares. It should not be used directly unless the
    user knows what he is doing.} 
  }
\section{Methods}{
\describe{

\item{\code{signature(model = "slinearModel")}}{
The method is specifically for system of linear models with the same
instruments and homoscedastic errors. It becomes SUR as a special case
when the instruments are the union of all regressors.
}
\item{\code{signature(model = "rslinearModel")}}{
This method is for restricted models that does not impose cross-equation
restrictions. With such restrictions 3SLS is not possible as we can
no longer write the model as a system of equations.
}
}}
\keyword{methods}

