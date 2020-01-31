\name{modelResponse-methods}
\docType{methods}
\alias{modelResponse}
\alias{modelResponse-methods}
\alias{modelResponse,linearModel-method}
\alias{modelResponse,slinearModel-method}
\alias{modelResponse,rslinearModel-method}
\alias{modelResponse,rlinearModel-method}
\title{ ~~ Methods for Function \code{modelResponse} in Package \pkg{momentfit} ~~}
\description{
Return the response vector in models with and without restrictions
}
\section{Methods}{
\describe{

\item{\code{signature(object = "linearModel")}}{
For linear models without restrictions on the coefficients.
}

\item{\code{signature(object = "slinearModel")}}{
For system of linear models without restrictions on the coefficients.
}

\item{\code{signature(object = "rslinearModel")}}{
For system of linear models with restrictions on the coefficients.
}

\item{\code{signature(object = "rlinearModel")}}{
For linear models with restrictions on the coefficients.
}
}}
\examples{
data(simData)

## Unrestricted model
model1 <- momentModel(y~x1+x2+x3, ~x2+x3+z1+z2, data=simData)
y <- modelResponse(model1)

## Restrictions change the response
R <- c("x2=2","x3=3")
rmodel1 <- restModel(model1, R)
rmodel1
restY <- modelResponse(rmodel1)

}
\keyword{methods}
