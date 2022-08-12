\name{solveGmm-methods}
\docType{methods}
\alias{solveGmm}
\alias{solveGmm-methods}
\alias{solveGmm,allNLModel,momentWeights-method}
\alias{solveGmm,linearModel,momentWeights-method}
\alias{solveGmm,slinearModel,sysMomentWeights-method}
\alias{solveGmm,rslinearModel,sysMomentWeights-method}
\alias{solveGmm,snonlinearModel,sysMomentWeights-method}
\alias{solveGmm,sfunctionModel,sysMomentWeights-method}
\alias{solveGmm,rnonlinearModel,momentWeights-method}
\title{ ~~ Methods for Function \code{solveGmm} in Package \pkg{momentfit} ~~}
\description{
The main function to get the GMM solution for a given weighting matrix.
}
\usage{
\S4method{solveGmm}{linearModel,momentWeights}(object, wObj, theta0=NULL,
...)

\S4method{solveGmm}{allNLModel,momentWeights}(object, wObj, theta0=NULL,
algo=c("optim","nlminb"), ...)

\S4method{solveGmm}{rnonlinearModel,momentWeights}(object, wObj, theta0=NULL,
...)

\S4method{solveGmm}{slinearModel,sysMomentWeights}(object, wObj, theta0=NULL)

\S4method{solveGmm}{rslinearModel,sysMomentWeights}(object, wObj, theta0=NULL)

\S4method{solveGmm}{snonlinearModel,sysMomentWeights}(object, wObj,
theta0=NULL, \dots)

\S4method{solveGmm}{sfunctionModel,sysMomentWeights}(object, wObj,
theta0=NULL, \dots)
}
\arguments{
  \item{object}{A moment-based model}
  \item{theta0}{The vector of coefficients for the starting values used
  in \code{\link{optim}}. If NULL, the starting values in the object if
  used. For system of equations, it is a list of vectors.}
  \item{wObj}{An object of class \code{"momentWeights"} or
    \code{"sysMomentWeights"}.}
  \item{algo}{The numerical algorithm to minimize the objective function.}
  \item{\dots}{Arguments to pass to \code{\link{optim}}.}
  }
\section{Methods}{
\describe{

\item{\code{signature(object = "allNLMoment", wObj = "momentWeights")}}{
Method to solve either nonlinear regressions or models in which moments
are computed with a function. The objective is minimized using optim.
}

\item{\code{signature(object = "rnonlinearModel", wObj = "momentWeights")}}{
Method to solve restricted nonlinear models. It computes the analytical
solution.
}

\item{\code{signature(object = "linearModel", wObj = "momentWeights")}}{
Method to solve linear models. It computes the analytical solution.
}

\item{\code{signature(object = "slinearModel", wObj = "sysMomentWeights")}}{
Method to solve system of linear models. It computes the analytical solution.
}

\item{\code{signature(object = "rslinearModel", wObj = "sysMomentWeights")}}{
Method to solve system of linear models in which restrictions have been
imposed on the coefficients. It computes the analytical solution.
}

\item{\code{signature(object = "slinearModel", wObj = "sysMomentWeights")}}{
Method to solve system of nonlinear models. The solution is obtained
with optim using the analytical derivatives.
}

}}
\value{
  A list with the following:
  \item{theta}{The vector of solution}
  \item{convergence}{convergence code. 0 means normal convergence. For
    higher numbers, see \code{\link{optim}}}
  }
\examples{
data(simData)
theta <- c(beta0=1,beta1=2)
model1 <- momentModel(y~x1, ~z1+z2, data=simData)

## A manual two-step GMM
w0 <- evalWeights(model1, w="ident")
theta0 <- solveGmm(model1, w0)$theta
w <- evalWeights(model1, theta0)
theta1 <- solveGmm(model1, w)$theta 
}

\keyword{methods}

