\name{evalGel-methods}
\docType{methods}
\alias{evalGel}
\alias{evalGel-methods}
\alias{evalGel,momentModel-method}
\title{ ~~ Methods for Function \code{evalGel} in Package \pkg{modelfit} ~~}
\description{
Method to simply evaluate a GEL model at a fixed coefficient vector. It
creates a \code{"gelfit"} object using that fixed vector. 
}
\usage{

\S4method{evalGel}{momentModel}(model, theta, lambda=NULL,
                                gelType="EL", rhoFct=NULL,
                                lamSlv=NULL, lControl=list(), \dots)
}
\arguments{
  \item{model}{An object of class \code{"momentModel"}.}
  \item{theta}{A vector of coefficients at which the model is estimated}
  \item{lambda}{The Lagrange multiplier vector. If not provided, the
    optimal vector is obtained for the given theta}
  \item{gelType}{The type of GEL. It is either \code{"EL"},
    \code{"ET"},  \code{"EEL"},  \code{"HD"},  \code{"ETEL"} or
    \code{"ETHD"}.}
  \item{rhoFct}{An alternative objective function for GEL. This argument
    is only used if we want to fit the model with a different GEL
    method. see \code{\link{rhoFct}}.} 
  \item{lamSlv}{An alternative solver for the Lagrange multiplier. By
    default, either \code{\link{Wu_lam}}, \code{\link{EEL_lam}},
    \code{\link{REEL_lam}} or \code{\link{getLambda}} is
    used.}
  \item{lControl}{A list of controls for the Lagrange multiplier
    algorithm.}
  \item{\dots}{Other arguments to pass. Not used for the moment.}
  }
\section{Methods}{
\describe{

\item{\code{signature(model = "momentModel")}}{
}
}}
\examples{
data(simData)
theta <- c(beta0=1,beta1=2)

## A linear model with optimal lambda
model1 <- momentModel(y~x1, ~z1+z2, data=simData)
evalGel(model1, c(1,1))

## A nonlinear model with fixed lambda
g <- y~beta0+x1^beta1
h <- ~z1+z2
model2 <- momentModel(g, h, c(beta0=1, beta1=2), data=simData)
evalGel(model2, theta=c(beta1=2, beta0=0.5), lambda=c(.1,.2,.3), gelType="ET")

}
\keyword{methods}
