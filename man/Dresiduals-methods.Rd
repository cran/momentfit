\name{Dresiduals-methods}
\docType{methods}
\alias{Dresiduals}
\alias{Dresiduals-methods}
\alias{Dresiduals,linearModel-method}
\alias{Dresiduals,nonlinearModel-method}
\alias{Dresiduals,sysModel-method}
\alias{Dresiduals,rsnonlinearModel-method}
\title{ ~~ Methods for Function \code{Dresiduals} in Package \pkg{Gmm} ~~}
\description{
It returns the matrix of derivatives of the residuals with respect to
the coefficients.
}
\section{Methods}{
\describe{

\item{\code{signature(object = "linearModel")}}{
%%  ~~describe this method here~~
}

\item{\code{signature(object = "nonlinearModel")}}{
%%  ~~describe this method here~~
}

\item{\code{signature(object = "rsnonlinearModel")}}{
%%  ~~describe this method here~~
}

\item{\code{signature(object = "sysMomentModel")}}{
%%  ~~describe this method here~~
}

}}

\examples{
data(simData)

theta <- c(beta0=1,beta1=2)
model1 <- momentModel(y~x1, ~z1+z2, data=simData)

Dresiduals(model1, theta)[1:3,]
}

\keyword{methods}