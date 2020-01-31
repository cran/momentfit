\name{tsls-methods}
\docType{methods}
\alias{tsls}
\alias{tsls-methods}
\alias{tsls,linearModel-method}
\alias{tsls,slinearModel-method}
\title{ ~~ Methods for Function \code{tsls} in Package \pkg{momentfit} ~~}
\description{
It estimates a linear model using two-stage least squares.
}
\usage{
\S4method{tsls}{linearModel}(model)

\S4method{tsls}{slinearModel}(model)
}
\arguments{
  \item{model}{An object of class \code{linearModel} or
    \code{slinearModel}.}
  }
\section{Methods}{
\describe{

\item{\code{signature(model = "linearModel")}}{
}

\item{\code{signature(model = "slinearModel")}}{
  2SLS for equation by equation estimation of a system of equations.
}

}}
\examples{
data(simData)
theta <- c(beta0=1,beta1=2)
model1 <- momentModel(y~x1, ~z1+z2, data=simData)
res <- tsls(model1)
summary(res)

## Econometrics, Fumio Hayashi (2000)
## Empirical exercises (b) and (c)
data(Griliches)
Griliches$YEAR <- as.factor(Griliches$YEAR)
model1  <- momentModel(LW~S+IQ+EXPR+TENURE+RNS+SMSA+YEAR-1,
                   ~S+EXPR+TENURE+RNS+SMSA+YEAR+MED+KWW+MRT+AGE-1,
                   data=Griliches, vcov="MDS")
res <- tsls(model1)
summary(res)
}
\keyword{methods}

