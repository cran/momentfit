\name{summary-methods}
\docType{methods}
\alias{summary}
\alias{summary-methods}
\alias{summary,gmmfit-method}
\alias{summary,gelfit-method}
\alias{summary,sgmmfit-method}
\title{ ~~ Methods for Function \code{summary} in Package \pkg{base} ~~}
\description{
Compute several results from a moment based model fit.
}
\usage{
\S4method{summary}{gmmfit}(object, testStrength=TRUE, \dots)

\S4method{summary}{gelfit}(object, \dots)

\S4method{summary}{sgmmfit}(object, testStrength=TRUE, \dots)
}
\arguments{
  \item{object}{A fit object from the package (GMM and GEL are the only
    methods for now)}
  \item{testStrength}{Should the first stage F-statistics be computed?}
  \item{\dots}{Other arguments to pass to \code{\link{vcov-methods}}}
  }
\section{Methods}{
\describe{

\item{\code{signature(object = "gmmfit")}}{
}

\item{\code{signature(object = "gmmfit")}}{
}


\item{\code{signature(object = "sgmmfit")}}{
}
}}
\examples{
data(simData)
theta <- c(beta0=1,beta1=2)
model1 <- momentModel(y~x1, ~z1+z2, data=simData)

res <- gmmFit(model1)
summary(res)

## Fixed and True Weights matrix
## Consider the moment of a normal distribution:
## Using the first three non centered moments

g <- function(theta, x)
{
mu <- theta[1]
sig2 <- theta[2]
m1 <- x-mu
m2 <- x^2-mu^2-sig2
m3 <- x^3-mu^3-3*mu*sig2
cbind(m1,m2,m3)
}

dg <- function(theta, x)
{
mu <- theta[1]
sig2 <- theta[2]
G <- matrix(c(-1,-2*mu,-3*mu^2-3*sig2, 0, -1, -3*mu),3,2)
}

x <- simData$x3
model <- momentModel(g, x, c(mu=.1, sig2=1.5), vcov="iid")
res1 <- gmmFit(model)
summary(res1)
## Same results (that's because the moment vcov is centered by default)
W <- solve(var(cbind(x,x^2,x^3)))
res2 <- gmmFit(model, weights=W)
res2
## If is therefore more efficient in this case to do the following:
summary(res2, breadOnly=TRUE)

}
\keyword{methods}

