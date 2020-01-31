\name{momentStrength-methods}
\docType{methods}
\alias{momentStrength}
\alias{momentStrength-methods}
\alias{momentStrength,functionModel-method}
\alias{momentStrength,formulaModel-method}
\alias{momentStrength,linearModel-method}
\alias{momentStrength,rlinearModel-method}
\alias{momentStrength,nonlinearModel-method}
\title{ ~~ Methods for Function \code{momentStrength} in Package \pkg{momentfit} ~~}
\description{
It produces measures of the strength of the moment conditons. 
}
\usage{
\S4method{momentStrength}{linearModel}(object, theta,
vcovType=c("OLS","HC","HAC","CL"))
}
\arguments{
  \item{object}{An object of class \code{"linearModel"}}
  \item{theta}{Coefficient vector at which the strength must be
    measured. It does not impact the measure for objects of class
    \code{linearModel}.}
  \item{vcovType}{Type of covariance matrix used to
    compute the F-test of the first-stage regression. For \code{HC},
    the function \code{\link{vcovHC}} is used with "HC1", and
    \code{\link{vcovHAC}} is used with the default setup is "HAC" is
    chosen. In \code{summary} for \code{gmmfit} objects, it is adjusted
    to the type of covariance that is set in the object. For type
  \code{CL}, clustered covariance matrix is computed. The options are
  the one included in the \code{vcovOptions} slot of the object (see
  \code{\link{meatCL}}). The object must have be defined with clusters
  for that to work. See \code{\link{momentModel}}.} 
}
\section{Methods}{
\describe{

\item{\code{signature(object = "functionModel")}}{
Not implemented yet. In that case, we want some measure of the rank of
the matrix of derivatives. 
}

\item{\code{signature(object = "formulaModel")}}{
Not implemented yet. In that case, we want some measure of the rank of
the matrix of derivatives. 
}

\item{\code{signature(object = "linearModel")}}{
It returns the F-test of the first stage regression. It is a measure of
the strength of the instruments.
}

\item{\code{signature(object = "rlinearModel")}}{
Returns nothing for now.
}

\item{\code{signature(object = "nonlinearModel")}}{
Not implemented yet. 
}
}}
\examples{
data(simData)

theta <- c(beta0=1,beta1=2)
model1 <- momentModel(y~x1, ~z1+z2, data=simData)
momentStrength(model1)
}
\keyword{methods}

