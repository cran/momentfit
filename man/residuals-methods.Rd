\name{residuals-methods}
\docType{methods}
\alias{residuals-methods}
\alias{residuals,ANY-method}
\alias{residuals,linearModel-method}
\alias{residuals,nonlinearModel-method}
\alias{residuals,sysModel-method}
\alias{residuals,rsysModel-method}
\alias{residuals,gmmfit-method}
\alias{residuals,gelfit-method}
\alias{residuals,sgmmfit-method}
\title{ ~~ Methods for Function \code{residuals} in Package \pkg{stats} ~~}
\description{
It computes the residual for a given coefficient vector, when the model
is a linear of nonlinear regression with instruments. The method can be
called on a \code{momentModel} class for a given coefficient theta or on a
\code{gmmfit} object.
}
\section{Methods}{
\describe{
\item{\code{signature(object = "rsysModel")}}{
}

\item{\code{signature(object = "linearModel")}}{
}

\item{\code{signature(object = "nonlinearModel")}}{
}

\item{\code{signature(object = "gmmfit")}}{
}
\item{\code{signature(object = "gelfit")}}{
}
\item{\code{signature(object = "sgmmfit")}}{
}
\item{\code{signature(object = "sysModel")}}{
}
}}
\examples{

x <- rchisq(200,5)
z1 <- rnorm(200)
z2 <- .2*x+rnorm(200)
y <- x+rnorm(200)
dat <- data.frame(y=y,z1=z1,x=x,z2=z2)
theta <- c(beta0=1,beta1=2)
model1 <- momentModel(y~x, ~z1+z2, data=dat)

## residuals for a given theta
e <- residuals(model1, theta)

## residuals of the fit
res <- gmmFit(model1)
e <- residuals(res)

}

\keyword{methods}
