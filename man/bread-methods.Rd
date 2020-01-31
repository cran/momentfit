\name{bread-methods}
\docType{methods}
\alias{bread}
\alias{bread-methods}
\alias{bread,gmmfit-method}
\alias{bread,sgmmfit-method}
\alias{bread,tsls-method}
\title{ ~~ Methods for Function \code{bread} in Package \pkg{sandwich} ~~}
\description{
It computes the bread in the sandwich representation of the covariance
matrix of the GMM estimator.  
}
\usage{
\S4method{bread}{gmmfit}(x, \dots)

\S4method{bread}{sgmmfit}(x, \dots)

\S4method{bread}{tsls}(x, \dots)
}
\arguments{
  \item{x}{GMM fit object}
  \item{\dots}{Arguments to pass to other methods}
  }
\section{Methods}{
\describe{

\item{\code{signature(x = "gmmfit")}}{
}

\item{\code{signature(x = "sgmmfit")}}{
}

\item{\code{signature(x = "tsls")}}{
}
}}
\examples{
data(simData)
theta <- c(beta0=1,beta1=2)
model1 <- momentModel(y~x1, ~z1+z2, data=simData)

res <- gmmFit(model1)
m <- meatGmm(res)
b <- bread(res)

## Sandwich vcov
b%*%m%*%b

## TSLS
model2 <- momentModel(y~x1, ~z1+z2, data=simData, vcov="iid")
res <- tsls(model2)
bread(res)


}
\keyword{methods}

