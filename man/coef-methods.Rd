\name{coef-methods}
\docType{methods}
\alias{coef-methods}
\alias{coef,gmmfit-method}
\alias{coef,gelfit-method}
\alias{coef,sgmmfit-method}
\alias{coef,rlinearModel-method}
\alias{coef,rslinearModel-method}
\alias{coef,rsnonlinearModel-method}
\alias{coef,rfunctionModel-method}
\alias{coef,rformulaModel-method}
\alias{coef,momentModel-method}
\alias{coef,rnonlinearModel-method}
\title{ ~~ Methods for Function \code{coef} in Package \pkg{stats} ~~}
\description{
  It extract the coefficient estimates of some moment-based models.
}
\section{Methods}{
\describe{
 
\item{\code{signature(object = "gmmfit")}}{
}

\item{\code{signature(object = "gelfit")}}{
}

\item{\code{signature(object = "sgmmfit")}}{
}

\item{\code{signature(object = "momentModel")}}{
}


\item{\code{signature(object = "rlinearModel")}}{
It gives the unrestricted representation of a restricted model. See examples.  
}

\item{\code{signature(object = "rslinearModel")}}{
It gives the unrestricted representation of a restricted model.  
}

\item{\code{signature(object = "rsnonlinearModel")}}{
It gives the unrestricted representation of a restricted model.  
}

\item{\code{signature(object = "rfunctionModel")}}{
It gives the unrestricted representation of a restricted model. See examples.  
}

\item{\code{signature(object = "rformulaModel")}}{
It gives the unrestricted representation of a restricted model. See examples.  
}

\item{\code{signature(object = "rnonlinearModel")}}{
It gives the unrestricted representation of a restricted nonlinear
model. 
}
}}
\examples{
data(simData)
model1 <- momentModel(y~x1+x2+x3+z1, ~x1+x2+z1+z2+z3+z4, data=simData)
res1 <- gmmFit(model1)
coef(res1)

### Restricted models
rmodel1 <- restModel(model1, R=c("x1=1", "x2=2*x3"))
res2 <- gmmFit(rmodel1)
res2
coef(rmodel1, coef(res2))
}
\keyword{methods}
