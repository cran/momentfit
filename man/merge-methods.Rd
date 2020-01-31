\name{merge-methods}
\docType{methods}
\alias{merge}
\alias{merge-methods}
\alias{merge,ANY,ANY-method}
\alias{merge,linearModel,linearModel-method}
\alias{merge,nonlinearModel,nonlinearModel-method}
\alias{merge,slinearModel,linearModel-method}
\alias{merge,snonlinearModel,nonlinearModel-method}
\title{ ~~ Methods for Function \code{merge} in Package \pkg{base} ~~}
\description{
It allows to merge \code{momentModel} classes into system objects. 
}
\usage{
\S4method{merge}{linearModel,linearModel}(x, y, \dots)

\S4method{merge}{nonlinearModel,nonlinearModel}(x, y, \dots)

\S4method{merge}{slinearModel,linearModel}(x, y, \dots)

\S4method{merge}{snonlinearModel,nonlinearModel}(x, y, \dots)

}
\arguments{
  \item{x}{An object on which the other objects are merged to.}
  \item{y}{An object to be merged to \code{x}.}
  \item{\dots}{Other objects of the same class as \code{y} to be merged
    to \code{x}.}
  }
\section{Methods}{
\describe{

\item{\code{signature(x = "linearModel", y = "linearModel")}}{
   Merging linear models into a system of equations.
}

\item{\code{signature(x = "nonlinearModel", y = "nonlinearModel")}}{
   Merging nonlinear models into a system of equations.
}

\item{\code{signature(x = "slinearModel", y = "linearModel")}}{
  Adding linear equations to a system of linear equations.
}

\item{\code{signature(x = "snonlinearModel", y = "nonlinearModel")}}{
Adding nonlinear equations to a system of nonlinear equations.
}
}}
\examples{
data(simData)
g1 <- y1~x1+x4; h1 <- ~z1+z2+z3+z4+x4
g2 <- y2~x1+x2+x3; h2 <- ~z1+z2+z3+z4+x3
g3 <- y3~x2+x3+x4; h3 <- ~z2+z3+z4+x3+x4
## Linear models
m1 <- momentModel(g1, h1, data=simData)
m2 <- momentModel(g2, h2, data=simData)
m3 <- momentModel(g3, h3, data=simData)

##
(sys1 <- merge(m1, m2))

## add an equation to the model

(sys2 <- merge(sys1, m3))

## want to get back the first?

sys2[1:2]

## Nonlinear (not really, just written as nonlinear)

nlg <- list(y1~theta0+theta1*x1+theta2*x4,
            y2~alpha0+alpha1*x1+alpha2*x2+alpha3*x3,
            y3~beta0+beta1*x2+beta2*x3+beta3*x4)
theta0 <- list(c(theta0=1,theta1=2,theta2=3),
              c(alpha0=1,alpha1=2,alpha2=3, alpha3=4),
              c(beta0=1,beta1=2,beta2=3,beta3=4))

nm1 <- momentModel(nlg[[1]], h1, theta0[[1]], data=simData)
nm2 <- momentModel(nlg[[2]], h2, theta0[[2]], data=simData)
nm3 <- momentModel(nlg[[3]], h3, theta0[[3]], data=simData)

merge(nm1, nm2, nm3)


}
\keyword{methods}

