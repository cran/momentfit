\name{momFct-methods}
\docType{methods}

\alias{momFct-methods}
\alias{momFct}
\alias{momFct,numeric,gelfit-method}

\title{Methods for Function \code{momFct} in Package \pkg{momentfit}}
\description{
 The methods computes the moment matrix. It is use to create special
 moment functions
}

\usage{
\S4method{momFct}{numeric,gelfit}(eta, object)
}

\arguments{
  \item{eta}{A vector that includes the coefficient and the Lagrange multipliers}
  \item{object}{An object of class \code{"gmmfit"}}
}

\section{Methods}{
\describe{

  \item{\code{signature(eta = "numeric", object = "gelfit")}}{
}
}}
\keyword{methods}
\keyword{misspecified}
