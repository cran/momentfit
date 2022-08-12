\name{modelDims-methods}
\docType{methods}
\alias{modelDims}
\alias{modelDims-methods}
\alias{modelDims,linearModel-method}
\alias{modelDims,rslinearModel-method}
\alias{modelDims,rsnonlinearModel-method}
\alias{modelDims,sysMomentModel-method}
\alias{modelDims,slinearModel-method}
\alias{modelDims,snonlinearModel-method}
\alias{modelDims,sfunctionModel-method}
\alias{modelDims,nonlinearModel-method}
\alias{modelDims,functionModel-method}
\alias{modelDims,formulaModel-method}
\alias{modelDims,rformulaModel-method}
\alias{modelDims,rlinearModel-method}
\alias{modelDims,rfunctionModel-method}
\alias{modelDims,rnonlinearModel-method}
\title{Methods for Function \code{modelDims}}
\description{
It extracts important information from the model. It is mostly used by
other methods when a \code{modelModel} has been modifed. An example is
when restrictions have been imposed on coefficients.
}
\section{Methods}{
\describe{

\item{\code{signature(object = "rlinearModel")}}{
}

\item{\code{signature(object = "rnonlinearModel")}}{
}

\item{\code{signature(object = "rfunctionModel")}}{
}

\item{\code{signature(object = "linearModel")}}{
}

\item{\code{signature(object = "nonlinearModel")}}{
}

\item{\code{signature(object = "functionModel")}}{
}

\item{\code{signature(object = "formulaModel")}}{
}

\item{\code{signature(object = "rformulaModel")}}{
}

\item{\code{signature(object = "slinearModel")}}{
}

\item{\code{signature(object = "rslinearModel")}}{
}

\item{\code{signature(object = "rsnonlinearModel")}}{
}

\item{\code{signature(object = "snonlinearModel")}}{
}

\item{\code{signature(object = "sfunctionModel")}}{
}
}}
\examples{
data(simData)

model1 <- momentModel(y~x1+x2, ~x2+z1+z2, data=simData)
modelDims(model1)

## Unrestricted model

rmodel1 <- restModel(model1, R=c("x1+x2=4"))
modelDims(rmodel1)


}

\keyword{methods}
