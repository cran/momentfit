\name{evalMoment-methods}
\docType{methods}
\alias{evalMoment}
\alias{evalMoment-methods}
\alias{evalMoment,functionModel-method}
\alias{evalMoment,formulaModel-method}
\alias{evalMoment,regModel-method}
\alias{evalMoment,sysModel-method}
\alias{evalMoment,rsysModel-method}
\title{ ~~ Methods for Function \code{evalMoment} in Package \pkg{momentfit} ~~}
\description{
Method to evaluate the moment matrix at a given coefficient vector.
}
\section{Methods}{
\describe{

\item{\code{signature(object = "functionModel")}}{
}

\item{\code{signature(object = "formulaModel")}}{
}

\item{\code{signature(object = "regModel")}}{
}

\item{\code{signature(object = "sysModel")}}{
}

\item{\code{signature(object = "rsysModel")}}{
}
}}
\examples{
data(simData)
theta <- c(1,1)
model1 <- momentModel(y~x1, ~z1+z2, data=simData)
gt <- evalMoment(model1, theta)

## A nonlinearGmm
g <- y~beta0+x1^beta1
h <- ~z1+z2
model2 <- momentModel(g, h, c(beta0=1, beta1=2), data=simData)
gt <- evalMoment(model2, c(beta0=1, beta1=2))

## A functionGmm
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
gt <- evalMoment(model3, c(beta0=1, beta1=2))
}
\keyword{methods}
