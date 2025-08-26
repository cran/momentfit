\name{quadra-methods}
\docType{methods}
\alias{quadra}
\alias{quadra-methods}
\alias{quadra,momentWeights,matrixORnumeric,matrixORnumeric-method}
\alias{quadra,momentWeights,matrixORnumeric,missing-method}
\alias{quadra,momentWeights,missing,missing-method}
\alias{quadra,sysMomentWeights,matrixORnumeric,matrixORnumeric-method}
\alias{quadra,sysMomentWeights,matrixORnumeric,missing-method}
\alias{quadra,sysMomentWeights,missing,missing-method}

\title{ ~~ Methods for Function \code{quadra} in Package \pkg{momentfit} ~~}
\description{
 ~~ Computes the quadratic form, where the center matrix is a class
 \code{momentWeights} object ~~
}
\usage{
\S4method{quadra}{momentWeights,missing,missing}(w, x, y, genInv=FALSE)

\S4method{quadra}{momentWeights,matrixORnumeric,missing}(w, x, y,
genInv=FALSE)

\S4method{quadra}{momentWeights,matrixORnumeric,matrixORnumeric}(w, x,
y, genInv=FALSE)


\S4method{quadra}{sysMomentWeights,matrixORnumeric,matrixORnumeric}(w,
x, y)

\S4method{quadra}{sysMomentWeights,matrixORnumeric,missing}(w, x, y)

\S4method{quadra}{sysMomentWeights,missing,missing}(w, x, y)


}
\arguments{
  \item{w}{An object of class \code{"momentWeights"}}
  \item{x}{A matrix or numeric vector}
  \item{y}{A matrix or numeric vector}
  \item{genInv}{Should we invert the center matrix using a generalized inverse?}
}

\section{Methods}{
\describe{
\item{\code{signature(w = "momentWeights", x = "matrixORnumeric",  y =
    "matrixORnumeric")}}{
It computes \eqn{x'Wy}, where \eqn{W} is the weighting matrix.
}
\item{\code{signature(w = "momentWeights", x = "matrixORnumeric",  y =
    "missing")}}{
It computes \eqn{x'Wx}, where \eqn{W} is the weighting matrix.  
}
\item{\code{signature(w = "momentWeights", x = "missing",  y =
    "missing")}}{
It computes \eqn{W}, where \eqn{W} is the weighting matrix.  When
\eqn{W} is the inverse of the covariance matrix of the moment
conditions, it is saved as either a QR decompisition, a Cholesky
decomposition or a covariance matrix into the \code{\link{momentWeights}}
object. The \code{quadra} method with no \code{y} and \code{x} is
therefore a way to invert it. The same applies to system of equations

}}}

\value{
  It returns a single numeric value.
}

\references{
  Courrieu P (2005), Fast Computation of Moore-Penrose Inverse Matrices.
  \emph{Neural Information Processing - Letters and Reviews}, \bold{8}(2), 25--29.
}

\examples{
data(simData)

theta <- c(beta0=1,beta1=2)
model1 <- momentModel(y~x1, ~z1+z2, data=simData)

gbar <- evalMoment(model1, theta)
gbar <- colMeans(gbar)

### Objective function of GMM with identity matrix
wObj <- evalWeights(model1, w="ident")
quadra(wObj, gbar)

### Objective function of GMM with efficient weights
wObj <- evalWeights(model1, theta)
quadra(wObj, gbar)

### Linearly dependent instruments

simData$z3 <- simData$z1+simData$z2
model2 <- momentModel(y~x1, ~z1+z2+z3, data=simData)
gbar2 <- evalMoment(model2, theta)
gbar2 <- colMeans(gbar2)

## A warning is printed about the singularity of the weighting matrix
wObj <- evalWeights(model2, theta)

## The regular inverse using the QR decomposition:
quadra(wObj)

## The regular inverse using the generalized inverse:
quadra(wObj, genInv=TRUE)


}
\keyword{methods}
