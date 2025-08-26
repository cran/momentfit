\name{minFit-methods}
\docType{methods}
\alias{minFit}
\alias{minFit-methods}
\alias{minFit,minAlgoNlm-method}
\alias{minFit,minAlgoStd-method}
\title{ Methods for Function \code{minFit} in Package \pkg{momentfit} ~~}
\description{
This is a method to minimize a function using an algorithm defined by
the \code{\link{minAlgo}} class. This is way of homogenizing the call of
minimization functions. It is used by \code{\link{solveGmm}} to offer
the possibility of using external solvers.
}
\usage{
\S4method{minFit}{minAlgoNlm}(object, start, fct, gr,
\dots)

\S4method{minFit}{minAlgoStd}(object, start, fct, gr,
\dots)
}

\arguments{
  \item{object}{A object of class \code{\link{minAlgo}} created by
  \code{\link{algoObj}}.} 
  \item{start}{A vector of starting values.}
  \item{fct}{The function to minimize.}    
  \item{gr}{An optional function that returns the gradient. The
    arguments of \code{fct} and \code{gr} must be identical.}
  \item{\dots}{Arguments to pass  the optimization algorithm and to the
  functions \code{fct} and \code{gr}.} 
}
\section{Methods}{
\describe{
\item{\code{signature(model = "minAlgoStd")}}{
This class includes all solvers that have a standard set of
  arguments. These arguments are the function, the gradient and the
  starting values (e.g. \code{\link{optim}}).
}

\item{\code{signature(model = "minAlgoNlm")}}{
This method is for solvers like \code{\link{nlm}}. The solver is quite
different from any other solvers, because the gradient is returned by
the main function as an attribute. That's why it needs a different
method to use it.
}
}}

\value{
  A list with the following elements:

  \item{solution}{The vector of solution for the minimization problem.}
  \item{value}{The value of the function at the solution.}
  \item{message}{The convergence message from the solver.}
  \item{convergence}{The convergence code from the solver.}

}

\examples{
f <- function(x, a=2, b=4) (x[1]-a)^2+(x[2]-b)^2
g <- function(x, a=2, b=4) c(2*(x[1]-a), 2*(x[2]-b))

## Using optim

algo1 <- algoObj("optim")
minFit(algo1, start=c(1,1), fct=f, gr=g, method="BFGS", b=6)

## Using nlm: very different structure, but same call

algo2 <- algoObj("nlm")
minFit(algo2, start=c(1,1), fct=f, gr=g, b=6)

}

\keyword{optimization}
