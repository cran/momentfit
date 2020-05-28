\name{evalGmm-methods}
\docType{methods}
\alias{evalGmm}
\alias{evalGmm-methods}
\alias{evalGmm,momentModel-method}
\alias{evalGmm,sysModel-method}
\title{ ~~ Methods for Function \code{evalGmm} in Package \pkg{modelfit} ~~}
\description{
Method to simply evaluate a GMM model at a fixed coefficient vector. It
creates a \code{"gmmfit"} object using that fixed vector. 
}
\usage{
\S4method{evalGmm}{momentModel}(model, theta, wObj=NULL, \dots)
\S4method{evalGmm}{sysModel}(model, theta, wObj=NULL, \dots)
}
\arguments{
  \item{model}{An object of class \code{"momentModel"}.}
  \item{theta}{A vector of coefficients at which the model is estimated}
  \item{wObj}{An object of class \code{"momentWeights"}. If not provided,
    the optimal weights based on the specification of the model
    evaluated at \code{theta} will be computed. }
  \item{\dots}{Other arguments to pass. Not used for the moment.}
  }
\section{Methods}{
\describe{

\item{\code{signature(model = "momentModel")}}{
}

\item{\code{signature(model = "sysModel")}}{
}
}}
\examples{
data(simData)
theta <- c(beta0=1,beta1=2)

## A linear model
model1 <- momentModel(y~x1, ~z1+z2, data=simData)
evalGmm(model1, c(1,1))

## A nonlinear model
g <- y~beta0+x1^beta1
h <- ~z1+z2
model2 <- momentModel(g, h, c(beta0=1, beta1=2), data=simData)
evalGmm(model2, theta=c(beta1=2, beta0=0.5))

## A function model
fct <- function(tet, x)
    {
        m1 <- (tet[1] - x)
        m2 <- (tet[2]^2 - (x - tet[1])^2)
        m3 <- x^3 - tet[1]*(tet[1]^2 + 3*tet[2]^2)
        f <- cbind(m1, m2, m3)
        return(f)
    }
dfct <- function(tet, x)
        {
        jacobian <- matrix(c( 1, 2*(-tet[1]+mean(x)), -3*tet[1]^2-3*tet[2]^2,0, 2*tet[2],
			   -6*tet[1]*tet[2]), nrow=3,ncol=2)
        return(jacobian)
        }
model3 <- momentModel(fct, simData$x3, theta0=c(beta0=1, beta1=2), grad=dfct)
evalGmm(model3, theta=c(beta1=.1, beta0=0.3))
}
\keyword{methods}
