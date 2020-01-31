\name{meatGmm-methods}
\docType{methods}
\alias{meatGmm}
\alias{meatGmm-methods}
\alias{meatGmm,gmmfit-method}
\alias{meatGmm,sgmmfit-method}
\alias{meatGmm,tsls-method}
\title{ ~~ Methods for Function \code{meatGmm} in Package \pkg{momentfit} ~~}
\description{
It computes the meat in the sandwich representation of the covariance
matrix of the GMM estimator.
}
\usage{
\S4method{meatGmm}{gmmfit}(object, robust=FALSE)

\S4method{meatGmm}{sgmmfit}(object, robust=FALSE)

\S4method{meatGmm}{tsls}(object, robust=FALSE)
}
\arguments{
  \item{object}{GMM fit object}
  \item{robust}{If \code{TRUE}, the meat is robust to the failure of the
    assumption that the weighting matrix is the inverse of the
    covariance matrix of the moment conditions. (see details)}
  }
\section{Methods}{
\describe{

  \item{\code{signature(object = "gmmfit")}}{
    General GMM fit.
}

\item{\code{signature(object = "tsls")}}{
  For model estimated by two-stage least squares.
}

\item{\code{signature(object = "sgmmfit")}}{
  For system of equations.
}
}}
\details{
If \code{robust=FALSE}, then the meat is \eqn{G'V^{-1}G}, where
\eqn{G} and \eqn{V} are respectively the sample mean of the derivatives
and the covariance matrix of the moment conditions. If it is
\code{TRUE}, the meat is \eqn{G'WVWG}, where \eqn{W} is the weighting
matrix.

For \code{tsls} objects, the function makes use of the QR representation
of the weighting matrix. It is simply possible to get the meat in a more
stable way. In that case, \eqn{W=(\sigma^2Z'Z/n)^{-1}}. If \code{robust}
is FALSE, \code{V} is assumed to be \eqn{\sigma^2Z'Z/n} which is the
inverse of the \code{bread}. Therefore, a sandwich covariance matrix
with \code{robust=FALSE} will result in a non-sandwich matrix.

For \code{sgmmfit}, the covariance is for the vectorized coefficient
vector of all equations. 
}
\examples{
data(simData)
theta <- c(beta0=1,beta1=2)
model1 <- momentModel(y~x1, ~z1+z2, data=simData)

res <- gmmFit(model1)
meatGmm(res)

## It is a slightly different because the weighting matrix
## is computed using the first step estimate and the covariance
## matrix of the moment conditions is based on the final estimate.
## They should, however, be asymptotically equivalent.

meatGmm(res, robust=TRUE)

## TSLS
res2 <- tsls(model1)

## Robust meat
meatGmm(res2, TRUE)

## It makes no difference is the model is assumed iid
model2 <- momentModel(y~x1, ~z1+z2, data=simData, vcov="iid")
res2 <- tsls(model2)
meatGmm(res2, FALSE)
meatGmm(res2, TRUE)


}
\keyword{methods}

