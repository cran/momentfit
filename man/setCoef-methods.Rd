\name{setCoef-methods}
\docType{methods}
\alias{setCoef}
\alias{setCoef-methods}
\alias{setCoef,momentModel-method}
\alias{setCoef,sysModel-method}

\title{Methods for Function \code{setCoef} in Package \pkg{momentfit} ~~}
\description{
  The method validates the coefficient theta and returns a coefficient
  object in a format that satisfies the moment model. 
}

\usage{
\S4method{setCoef}{momentModel}(model, theta)
\S4method{setCoef}{sysModel}(model, theta)

}
\arguments{
  \item{model}{A moment model object.}  
  \item{theta}{A coefficient object. The type depends on the model
    object. See the examples below.} 
}

\section{Methods}{
\describe{
  \item{\code{signature(object = "momentModel")}}{
    Methods for all single equation models including the restricted
    ones.    
  }

  \item{\code{signature(object = "sysModel")}}{
    Methods for all system of equations models including the restricted
    ones.  
  }
}}
\examples{
### A few system of equation models:
data(simData)
h <- list(~z1+z2+z3, ~x3+z1+z2+z3+z4, ~x3+x4+z1+z2+z3)
nlg <- list(Supply=y1~theta0+theta1*x1+theta2*z2,
            Demand1=y2~alpha0+alpha1*x1+alpha2*x2+alpha3*x3,
            Demand2=y3~beta0+beta1*x3+beta2*x4+beta3*z1)
g <- list(Supply=y1~x1+z2, Demand1=y2~x1+x2+x3, Demand2=y3~x3+x4+z1)
theta0 <- list(c(theta0=1,theta1=2,theta2=3),
               c(alpha0=1,alpha1=2,alpha2=3, alpha3=4),
               c(beta0=1,beta1=2,beta2=3,beta3=4))
nlin <- sysMomentModel(nlg, h, theta0, data=simData)
lin <- sysMomentModel(g, h, data=simData)

### from numeric vector to the proper format with names:
setCoef(nlin, 1:11)

### reorder the equation and name the coefficients
setCoef(nlin, list(Demand1=1:4, Supply=1:3, Demand2=1:4))

### reorder the coefficint to match the order in the model
tet <- do.call("c", theta0)
set.seed(112233)
setCoef(nlin, tet[sample(11)])

### It validates length and names and provide source of errors
\dontrun{
setCoef(nlin, list(Demand1=1:4, Supply=1:2, Demand2=1:4))
names(tet)[4] <- "gamma3"
setCoef(nlin, tet)
setCoef(nlin, list(Demand1=1:4, Supply=1:3, Demand4=1:4))
}

### a single equation model
single <- momentModel(nlg[[1]], h[[1]], theta0[[1]], data=simData)
setCoef(single, c(theta1=4, theta0=6, theta2=8))
setCoef(single, 1:3)
}
\keyword{methods}
