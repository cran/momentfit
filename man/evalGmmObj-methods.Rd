\name{evalGmmObj-methods}
\docType{methods}
\alias{evalGmmObj}
\alias{evalGmmObj-methods}
\alias{evalGmmObj,momentModel,numeric,momentWeights-method}
\alias{evalGmmObj,sysModel,list,sysMomentWeights-method}
\title{ ~~ Methods for Function \code{evalGmmObj} in Package \pkg{momentfit} ~~}
\description{
 ~~ Methods to compute the GMM objective function. ~~
}
\usage{
\S4method{evalGmmObj}{momentModel,numeric,momentWeights}(object, theta,
wObj, \dots)

\S4method{evalGmmObj}{sysModel,list,sysMomentWeights}(object, theta,
wObj, \dots)
}
\arguments{
  \item{object}{An object of class \code{"momentModel"},
    or \code{"sysMomentModels"}.}
  \item{theta}{The vector for coefficients for single equation, or a
    list of vector for system of equations.}
  \item{wObj}{An object of class \code{"momentWeights"} or
    \code{"sysMomentWeights"}.}
  \item{\dots}{Arguments to pass to other methods}
  }
\section{Methods}{
\describe{

\item{\code{signature(object = "momentModel", theta = "numeric",
    wObj = "momentWeights")}}{
}

\item{\code{signature(object = "sysModel", theta = "list",
    wObj = "sysMomentWeights")}}{
}
}}
\examples{
data(simData)

theta <- c(beta0=1,beta1=2)
model1 <- momentModel(y~x1, ~z1+z2, data=simData)
w <- evalWeights(model1, theta)
evalGmmObj(model1, theta, w)

}

\keyword{methods}
