\name{algoObj}

\alias{algoObj}
	
\title{Constructor for \code{\link{minAlgo}} classes}

\description{
This function creates an object that defines minimization solvers. The
purpose is to homogenize the call of optimization solvers. These objects
can be used by the method \code{\link{minFit}} to call the associated
solvers using the same arguments, and to return the solution using
a list with the same names.
}

\usage{
algoObj(algo, start, fct, grad, solution, value, message, convergence)  
}

\arguments{
\item{algo}{The name of the solver function to be called in character
  format. All arguments for the solver \code{\link{optim}},
  \code{\link{nlminb}}, \code{\link{constrOptim}} and \code{\link{nlm}}
  are determined automatically. The other arguments are only needed for
  solvers coming from other packages.}

\item{start}{The name of the argument representing the starting value
  in character format.}  

\item{fct}{The name of the argument representing the function to
  minimize in character format.}  

\item{grad}{The name of the argument representing the gradient in
  character format.}   

\item{solution}{The name of the element of the list returned by the
  solver that represents the solution, in character format.}   

\item{value}{The name of the element of the list returned by the
  solver that represents the value of the function at the solution, in
  character format.}  

\item{message}{The name of the element of the list returned by the
  solver that represents the convergence message, in character format.}

\item{convergence}{The name of the element of the list returned by the
  solver that represents the convergence code, in character format.}
}

\value{
  An object of class \code{\link{minAlgo}}.
 }

 \seealso{
\code{\link{minFit}} for examples on how to use this class object.
 }
 
\examples{
## The optim algorithm:
algo1 <- algoObj("optim")

## The nlminb algorithm:
algo2 <- algoObj("nlminb")

## Defining the algorithm lbfgs from the nloptr package

\dontrun{
algo3 <- algoObj(algo="lbfgs", start="x0", fct="fn",
                 grad="gr", solution="par", value="value",
                 message="message", convergence="convergence")
}

}

