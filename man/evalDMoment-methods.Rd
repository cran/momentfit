\name{evalDMoment-methods}
\docType{methods}
\alias{evalDMoment}
\alias{evalDMoment-methods}
\alias{evalDMoment,functionModel-method}
\alias{evalDMoment,rfunctionModel-method}
\alias{evalDMoment,formulaModel-method}
\alias{evalDMoment,rformulaModel-method}
\alias{evalDMoment,rnonlinearModel-method}
\alias{evalDMoment,sysModel-method}
\alias{evalDMoment,rslinearModel-method}
\alias{evalDMoment,rsnonlinearModel-method}
\alias{evalDMoment,regModel-method}
\title{ ~~ Methods for Function \code{evalDMoment} in Package \pkg{momentfit} ~~}
\description{
It computes the matrix of derivatives of the sample moments with respect
to the coefficients.
}

\usage{
\S4method{evalDMoment}{functionModel}(object, theta, impProb=NULL,
lambda=NULL)

\S4method{evalDMoment}{rfunctionModel}(object, theta, impProb=NULL,
lambda=NULL)

\S4method{evalDMoment}{rnonlinearModel}(object, theta, impProb=NULL,
lambda=NULL)

\S4method{evalDMoment}{formulaModel}(object, theta, impProb=NULL,
lambda=NULL)

\S4method{evalDMoment}{rformulaModel}(object, theta, impProb=NULL,
lambda=NULL)

\S4method{evalDMoment}{regModel}(object, theta, impProb=NULL,
lambda=NULL)

\S4method{evalDMoment}{sysModel}(object, theta)

\S4method{evalDMoment}{rslinearModel}(object, theta)

\S4method{evalDMoment}{rsnonlinearModel}(object, theta, impProb=NULL)
}

\arguments{
  \item{object}{An model object}
  \item{theta}{A numerical vector of coefficients} 
  \item{impProb}{If a vector of implied probablities is provided, the
    sample means are computed using them. If not provided, the means are
    computed using the uniform weight}
  \item{lambda}{A vector of Lagrange multipliers associated with the
  moment conditions. Its length must therefore match the number of
  conditions. See details below.}
}
  
\section{Methods}{
\describe{

\item{\code{signature(object = "functionModel")}}{
}

\item{\code{signature(object = "rfunctionModel")}}{
The theta vector must match the number of coefficients in the restricted
model.
}

\item{\code{signature(object = "formulaModel")}}{
}

\item{\code{signature(object = "rformulaModel")}}{
The theta vector must match the number of coefficients in the restricted
model.
}

\item{\code{signature(object = "regModel")}}{
}

\item{\code{signature(object = "sysModel")}}{
}

\item{\code{signature(object = "rslinearModel")}}{
}
}}

\details{
Without the argument \code{lambda}, the method returns a \eqn{q \times
  k} matrix, where \eqn{k} is the number of coefficients, and \eqn{q} is
the number of moment conditions. That matrix is the derivative of the
sample mean of the moments with respect to the coefficient. 

If \code{lambda} is provided, the method returns an \eqn{n \times k}
  matrix, where \eqn{n} is the sample size. The ith row is
  \eqn{G_i'\lambda}, where $G_i$ is the derivative of the moment
  function evaluated at the ith observation. For now, this option is
  used to compute robust-to-misspecified standard errors of GEL
  estimators.
}


\examples{
data(simData)
theta <- c(1,1)
model1 <- momentModel(y~x1, ~z1+z2, data=simData)
G <- evalDMoment(model1, theta)

## A nonlinearModel
g <- y~beta0+x1^beta1
h <- ~z1+z2
model2 <- momentModel(g, h, c(beta0=1, beta1=2), data=simData)
G <- evalDMoment(model2, c(beta0=1, beta1=2))

## A functionModel
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
X <- rnorm(200)
model3 <- momentModel(fct, X, theta0=c(beta0=1, beta1=2), grad=dfct)
G <- evalDMoment(model3, c(beta0=1, beta1=2))
}

\keyword{methods}
