\name{vcovHAC-methods}
\docType{methods}
\alias{vcovHAC-methods}
\alias{vcovHAC,momentModel-method}
\title{ ~~ Methods for Function \code{vcovHAC} in Package \pkg{sandwich} ~~}
\description{
 Methods to compute the HAC covariance matrix of the moment matrix ~~
}
\section{Methods}{
\describe{

\item{\code{signature(x = "momentModel")}}{
}
}}

\examples{
data(simData)
theta <- c(beta0=1,beta1=2)
model1 <- momentModel(y~x1, ~z1+z2, data=simData)

# a warning is given if the model is not set as being HAC
vcovHAC(model1, theta)

model1 <- momentModel(y~x1, ~z1+z2, data=simData, vcov="HAC",vcovOptions=list(kernel="B"))
vcovHAC(model1, theta)
}

\keyword{methods}

