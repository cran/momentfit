\name{evalGelObj-methods}
\docType{methods}
\alias{evalGelObj}
\alias{evalGelObj-methods}
\alias{evalGelObj,momentModel,numeric,numeric-method}
\title{ ~~ Methods for Function \code{evalGelObj} in Package \pkg{Gmm} ~~}
\description{
 ~~ Methods to compute the GEL objective function. ~~
}
\usage{
\S4method{evalGelObj}{momentModel,numeric,numeric}(object, theta,
                                                   lambda, gelType,
                                                   rhoFct=NULL, \dots)
}
\arguments{
  \item{object}{An object of class \code{"momentModel"}}
  \item{theta}{The vector for coefficients.}
  \item{lambda}{Vector of Lagrange multiplier.}
  \item{gelType}{The type of GEL. It is either \code{"EL"},
    \code{"ET"},  \code{"EEL"},  \code{"HD"},  \code{"ETEL"} or
    \code{"ETHD"}.}
  \item{rhoFct}{An alternative objective function for GEL. This argument
    is only used if we want to fit the model with a different GEL
    method. see \code{\link{rhoFct}}.} 
  \item{\dots}{Arguments to pass to other methods}
  }
\section{Methods}{
\describe{

\item{\code{signature(object = "momentModel", theta = "numeric",
    lambda = "numeric")}}{
}
}}
\examples{
data(simData)

theta <- c(beta0=1,beta1=2)
model1 <- momentModel(y~x1, ~z1+z2, data=simData)
evalGelObj(model1, theta, c(.2,.3,.4), gelType="EL")

}

\keyword{methods}
