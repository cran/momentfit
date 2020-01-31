\name{vcov-methods}
\docType{methods}
\alias{vcov-methods}
\alias{vcov,gmmfit-method}
\alias{vcov,momentModel-method}
\alias{vcov,sysModel-method}
\alias{vcov,gelfit-method}
\alias{vcov,sgmmfit-method}
\alias{vcov,tsls-method}
\title{ ~~ Methods for Function \code{vcov} in Package \pkg{stats} ~~}
\description{
Computes the covariance matrix of the coefficient estimated by GMM or GEL. 
}
\usage{
\S4method{vcov}{gmmfit}(object, sandwich=NULL, df.adj=FALSE,
breadOnly=FALSE, modelVcov=NULL)

\S4method{vcov}{sgmmfit}(object, sandwich=NULL, df.adj=FALSE,
breadOnly=FALSE, modelVcov=NULL)

\S4method{vcov}{tsls}(object, sandwich=TRUE, df.adj=FALSE)

\S4method{vcov}{gelfit}(object, withImpProb=FALSE, tol=1e-10,
                        robToMiss=FALSE)

\S4method{vcov}{momentModel}(object, theta)

\S4method{vcov}{sysModel}(object, theta)
}
\arguments{
  \item{object}{A fitted model or a model, For fitted models, it
    computes the covariance matrix of the estimators. For models, it
    computes the covariance matrix of the moment conditions, in which
    case, the coefficient vector must be provided. }
  \item{theta}{Coefficient vector to compute the covariance matrix of
    the moment conditions.}
  \item{sandwich}{Should we compute the sandwich covariance matrix. This is
    only necessary if the weighting matrix is not the optimal one, or if
    we think it is a bad estimate of it. If \code{NULL}, it will be set
    to "TRUE" for One-Step GMM, which includes just-identified GMM like
    IV, and "FALSE" otherwise. }
  \item{df.adj}{Should we adjust for degrees of freedom. If \code{TRUE}
    the covariance matrix is multiplied by \code{n/(n-k)}, where
    \code{n} is the sample size and \code{k} is the number of
    coefficients. For heteroscedastic robust covariance matrix,
    adjusting is equivalent to computing HC1 while not adjusting is
    HC0.}
  \item{breadOnly}{If \code{TRUE}, the covariance matrix is set to the
    bread (see details below).}
  \item{modelVcov}{Should be one of "iid", "MDS" or "HAC". It is meant
    to change the way the variance of the moments is computed. If it is
    set to a different specification included in the model,
    \code{sandwich} is set to \code{TRUE}. }
  \item{withImpProb}{Should we compute the moments with the implied
    probabilities}
  \item{tol}{Any diagonal less than \code{"tol"} is set to tol}
  \item{robToMiss}{Should we compute a covariance matrix that is robust
    to misspecification?}
}
\section{Methods}{
\describe{

\item{\code{signature(object = "gmmfit")}}{
For any model estimated by any GMM methods.  
}

\item{\code{signature(object = "gelfit")}}{
For any model estimated by any GMM methods.  
}

\item{\code{signature(object = "sgmmfit")}}{
For any system of equations estimated by any GMM methods.  
}
}}
\details{
If \code{sandwich=FALSE}, then it returns \eqn{(G'V^{-1}G)^{-1}/n}, where
\eqn{G} and \eqn{V} are respectively the matrix of average derivatives
and the covariance matrix of the moment conditions. If it is
\code{TRUE}, it returns \eqn{(G'WG)^{-1}G'WVWG(G'WG)^{-1}/n},
where \eqn{W} is the weighting matrix used to obtain the vector of
estimates.

If \code{breadOnly=TRUE}, it returns \eqn{(G'WG)^{-1}/n},
where the value of \eqn{W} depends on the type of GMM. For two-step GMM,
it is the first step weighting matrix, for one-step GMM, it is either
the identity matrix or the fixed weighting matrix that was provided when
\code{\link{gmmFit}} was called, for iterative GMM, it is the weighting
matrix used in the last step. For CUE, the result is identical to
\code{sandwich=FALSE} and \code{beadOnly=FALSE}, because the
weighting and coefficient estimates are obtained simultaneously, which
makes \eqn{W} identical to \eqn{V}.

\code{breadOnly=TRUE} should therefore be used with caution because it
will produce valid standard errors only if the weighting matrix
converges to the the inverse of the covariance matrix of the moment
conditions.

For \code{"tsls"} objects, \code{sandwich} is TRUE by default. If we
assume that the error term is iid, then setting it to FALSE to result in
the usual \eqn{\sigma^2(\hat{X}'\hat{X})^{-1}} covariance matrix. If
\code{FALSE}, it returns a robust covariance matrix determined by the
value of \code{vcov} in the \code{momentModel}. 
}
\examples{
data(simData)
theta <- c(beta0=1,beta1=2)
model1 <- momentModel(y~x1, ~z1+z2, data=simData)

## optimal matrix
res <- gmmFit(model1)
vcov(res)

## not the optimal matrix
res <- gmmFit(model1, weights=diag(3))
vcov(res, TRUE)

## Model with heteroscedasticity
## MDS is for models with no autocorrelation.
## No restrictions are imposed on the structure of the
## variance of the moment conditions
model2 <- momentModel(y~x1, ~z1+z2, data=simData, vcov="MDS")
res <- tsls(model2)

## HC0 type of robust variance
vcov(res, sandwich=TRUE)
## HC1 type of robust variance
vcov(res, sandwich=TRUE, df.adj=TRUE)

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
## the option breadOnly of summary() is passed to vcov()
summary(res2, breadOnly=TRUE)
}
\keyword{methods}

