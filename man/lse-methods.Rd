\name{lse-methods}
\docType{methods}
\alias{lse-methods}
\alias{lse}
\alias{lse,linearModel-method}
\title{Least Squares Methods for Moment Models}
\description{
It estimates models defined in the package by least squares. At the
moment, it only applies to \code{linearModel} objects and it is
estimated using \code{\link{lm}}.
}
\section{Methods}{
\describe{
\item{\code{signature(model = "linearModel")}}{
It ignores the instruments and simply fits the linear model with LSE.
}
}}

\examples{
data(simData)
mod <- momentModel(y~x1, ~z1+z2, vcov="MDS", data=simData)
lse(mod)
}

\keyword{methods}
\keyword{LSE}
