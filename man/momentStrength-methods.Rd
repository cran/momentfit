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
It produces measures of the strength of the moment conditions.
}
\usage{
\S4method{momentStrength}{linearModel}(object, theta) 
}
\arguments{
  \item{object}{An object of class \code{"linearModel"}}
  \item{theta}{Coefficient vector at which the strength must be
    measured. It does not impact the measure for objects of class
    \code{linearModel}.}
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

\details{
For now, the method only exists for linear models. It returns
the F-statistics from the first stage regression. The type of covariance
matrix used to compute the statistics depends on the specification of
the model. If the argument \code{vcov} of the model is set to
\code{"iid"}, a non robust estimator is used. If it is set to
\code{"MDS"}, \code{"HAC"}, or \code{"CL"}, the appropriate robust
estimator is used. To use a different type, use the method \code{update}
to change the argument \code{vcov} of the model object. See the
example below.
}

\examples{
data(simData)

theta <- c(beta0=1,beta1=2)
model1 <- momentModel(y~x1, ~z1+z2, data=simData, vcov="iid")
momentStrength(model1)
## changing the type of vcov to get robust tests
momentStrength(update(model1, vcov="MDS"))
}
\keyword{methods}

