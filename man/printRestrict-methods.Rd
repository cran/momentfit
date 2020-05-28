\name{printRestrict-methods}
\docType{methods}
\alias{printRestrict}
\alias{printRestrict-methods}
\alias{printRestrict,rlinearModel-method}
\alias{printRestrict,rformulaModel-method}
\alias{printRestrict,rslinearModel-method}
\alias{printRestrict,rnonlinearModel-method}
\alias{printRestrict,rsnonlinearModel-method}
\alias{printRestrict,rfunctionModel-method}
\title{ ~~ Methods for Function \code{printRestrict} in Package \pkg{momentfit} ~~}
\description{
It prints the detailed restrictions imposed on \code{"momentModel"} classes.
}
\section{Methods}{
\describe{

\item{\code{signature(object = "rgelModels")}}{
}  
  
\item{\code{signature(object = "rlinearModel")}}{
}

\item{\code{signature(object = "rnonlinearModel")}}{
}

\item{\code{signature(object = "rfunctionModel")}}{
}

\item{\code{signature(object = "rformulaModel")}}{
}

\item{\code{signature(object = "rslinearModel")}}{
}

\item{\code{signature(object = "rsnonlinearModel")}}{
}
}}
\examples{
data(simData)
theta <- c(beta0=1,beta1=2)

## Unrestricted model
model1 <- momentModel(y~x1+x2+x3+z1, ~x1+x2+z1+z2+z3+z4, data=simData)

## restricted model
R <- matrix(c(1,1,0,0,0,0,0,2,0,0,0,0,0,1,-1),3,5, byrow=TRUE)
q <- c(0,1,3)
rmodel1 <- restModel(model1, R, q)
printRestrict(rmodel1)
}
\keyword{methods}
