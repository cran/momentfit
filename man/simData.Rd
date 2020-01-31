\name{simData}
\alias{simData}
\docType{data}
\title{
Simulated data.
}
\description{
This dataset is used in several documentation files to illustrate the
different functionality of the package.
}
\usage{data("simData")}
\format{
  A data frame with 50 observations on the following 12 variables. See
  the examples for the method used to generate them.
  \describe{
    \item{\code{y}}{a numeric vector}
    \item{\code{y1}}{a numeric vector}
    \item{\code{y3}}{a numeric vector}
    \item{\code{y2}}{a numeric vector}
    \item{\code{z1}}{a numeric vector}
    \item{\code{x1}}{a numeric vector}
    \item{\code{z2}}{a numeric vector}
    \item{\code{x2}}{a numeric vector}
    \item{\code{z3}}{a numeric vector}
    \item{\code{x3}}{a numeric vector}
    \item{\code{x4}}{a numeric vector}
    \item{\code{z4}}{a numeric vector}
    \item{\code{z5}}{a numeric vector}
  }
}

\examples{
# Here is how the data was simulated
set.seed(1122)
n <- 50
x1 <- rchisq(n,5)
x2 <- rchisq(n,5)
x3 <- rnorm(n)
x4 <- rnorm(n)
z1 <- .2*x1+rnorm(n)
z2 <- .2*x2+rnorm(n)
z3 <- rnorm(n)
z4 <- rnorm(n)
z5 <- rnorm(n)
y <- y1 <- x1+rnorm(n)
y2 <- 2*x1+rnorm(n)
y3 <- 0.5*x2+rnorm(n)
simData <- data.frame(y=y, y1=y1,y3=y3,y2=y2, z1=z1,x1=x1,z2=z2,x2=x2,z3=z3,x3=x3,
                  x4=x4,z4=z4,z5=z5)
}
\keyword{datasets}
