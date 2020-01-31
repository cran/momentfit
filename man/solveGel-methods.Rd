\name{solveGel-methods}
\docType{methods}
\alias{solveGel}
\alias{solveGel-methods}
\alias{solveGel,momentModel-method}
\title{ ~~ Methods for Function \code{solveGel} in Package \pkg{momentfit} ~~}
\description{
It fits a moment-based model using GEL methods.
}
\usage{
\S4method{solveGel}{momentModel}(object, gelType="EL", theta0=NULL,
                                  lambda0=NULL, lamSlv=NULL,
                                  coefSlv=c("optim","nlminb","constrOptim"),
                                  rhoFct=NULL, 
                                  lControl=list(), tControl=list())
}
\arguments{
  \item{object}{An object of class \code{"gelModels"}}
  \item{gelType}{The type of GEL. It is either \code{"EL"},
    \code{"ET"},  \code{"EEL"},  \code{"HD"},  \code{"ETEL"} or
    \code{"ETHD"}.}
  \item{theta0}{The vector of coefficients for the starting values used
    in minimization algorithm. If NULL, the starting values in the object is
    used. For linear models, it must be provided because
    \code{"linearGel"} objects do not have a theta0 slot.}
  \item{lambda0}{The \eqn{q \times 1} vector of starting values for the
    Lagrange multipliers. By default a zero vector is used.}
  \item{lamSlv}{An alternative solver for the Lagrange multiplier. By
    default, either \code{\link{Wu_lam}}, \code{\link{EEL_lam}},
    \code{\link{REEL_lam}} or \code{\link{getLambda}} is
    used.}
  \item{coefSlv}{Minimization solver for the coefficient vector.}	    
  \item{rhoFct}{An alternative objective function for GEL. This argument
    is only used if we want to fit the model with a different GEL
    method. see \code{\link{rhoFct}}.} 
  \item{lControl}{A list of controls for the Lagrange multiplier
    algorithm.}
  \item{tControl}{A list of controls for the coefficient algorithm.}
  }
\section{Methods}{
\describe{

\item{\code{signature(object = "momentModel")}}{
The method applies to all GEL classes.
}

}}
\value{
  A list with the following:
  \item{theta}{The vector of solution}
  \item{lambda}{The vector of Lagrange multiplier}
  \item{lconvergence}{convergence code for the Lagrange multiplier. 0
    means normal convergence.}
  \item{convergence}{convergence code for the coefficients. 0 means
    normal convergence. For higher numbers, see \code{\link{optim}},
    \code{\link{constrOptim}} or  \code{\link{nlminb}}}
  }
\examples{
data(simData)
model1 <- momentModel(y~x1, ~z1+z2, data=simData)

## Get a good starting value
theta0 <- gmmFit(model1)@theta

## EL by default, with Wu algorithm
res2 <- solveGel(model1, theta0=theta0)

## Change solver parameters
res3 <- solveGel(model1, theta0=theta0,
                 tControl=list(method="Nelder", control=list(maxit=2000)))

}

\keyword{methods}

