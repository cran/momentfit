\name{gelFit-methods}
\docType{methods}
\alias{gelFit}
\alias{gelFit-methods}
\alias{gelFit,momentModel-method}
\alias{gelFit,rmomentModel-method}

\title{ ~~ Methods for Function \code{gelFit} in Package \pkg{momentfit} ~~}
\description{
Method to fit a model using GEL, from an object of class
\code{"momentModel"} or its restricted counterpart.
}
\usage{
\S4method{gelFit}{momentModel}(model, gelType="EL", rhoFct=NULL,
              initTheta=c("gmm", "modelTheta0"), theta0=NULL,
              lambda0=NULL, vcov=FALSE, \dots)

\S4method{gelFit}{rmomentModel}(model, gelType="EL", rhoFct=NULL,
              initTheta=c("gmm", "modelTheta0"), theta0=NULL,
              lambda0=NULL, vcov=FALSE, ...)
}
\arguments{
  \item{model}{A model class object}
  \item{gelType}{The type of GEL. It is either \code{"EL"},
              \code{"ET"},  \code{"EEL"},  \code{"HD"},  \code{"ETEL"} or
               \code{"ETHD"}.}
  \item{rhoFct}{An alternative objective function for GEL. This argument
    is only used if we want to fit the model with a different GEL
    method. see \code{\link{rhoFct}}.}
  \item{initTheta}{Method to obtain the starting values for the
    coefficient vector. By default the GMM estimate with identity matrix
    is used. The second argument means that the theta0 of the
    object, if any, should be used.}  
  \item{theta0}{An optional initial vector for \code{\link{optim}} when
    the model is nonlinear. If provided, the argument \code{"initTheta"}
    is ignored.}
  \item{lambda0}{Manual starting values for the Lagrange
    multiplier. By default, it is a vector of zeros.}
  \item{vcov}{Should the method computes the covariance matrices of the
    coefficients and Lagrange multipliers.}
  \item{\dots}{Arguments to pass to other methods (mostly the
    optimization algorithm)}  
  }
\section{Methods}{
\describe{
\item{\code{signature(model = "momentModel")}}{
The main method for all moment-based  models.
}

\item{\code{signature(model = "rmomentModel")}}{
The main method for all restricted moment-based  models.
}

}}

\examples{
data(simData)

theta <- c(beta0=1,beta1=2)
model1 <- momentModel(y~x1, ~z1+z2, data=simData)

## EL estimate
res1 <- gelFit(model1)
res1

## ET estimate
res2 <- gelFit(model1, gelType="ET")
res2

## Restricted models by EL
## using the Brent method
R <- matrix(c(0,1), ncol=2)
q <- 2
rmodel1 <- restModel(model1, R, q)
gelFit(rmodel1, tControl=list(method="Brent", lower=-10, upper=10))

}
\keyword{methods}
