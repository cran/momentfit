\name{model.matrix-methods}
\docType{methods}
\alias{model.matrix}
\alias{model.matrix-methods}
\alias{model.matrix,linearModel-method}
\alias{model.matrix,nonlinearModel-method}
\alias{model.matrix,rlinearModel-method}
\alias{model.matrix,slinearModel-method}
\alias{model.matrix,rslinearModel-method}
\alias{model.matrix,rsnonlinearModel-method}
\alias{model.matrix,snonlinearModel-method}
\title{ ~~ Methods for Function \code{model.matrix} in Package \pkg{stats} ~~}
\description{
 Model matrix form \code{momentModel}. It returns the matrix of regressors
 or the instruments. In restricted models, it returns the reduced matrix
 of regressors.
}
\usage{
\S4method{model.matrix}{linearModel}(object,
type=c("regressors","instruments"))
\S4method{model.matrix}{rlinearModel}(object,
type=c("regressors","instruments"))
\S4method{model.matrix}{nonlinearModel}(object,
type=c("regressors","instruments"))
\S4method{model.matrix}{slinearModel}(object,
type=c("regressors","instruments"))
\S4method{model.matrix}{rslinearModel}(object,
type=c("regressors","instruments"))
\S4method{model.matrix}{rsnonlinearModel}(object,
type=c("regressors","instruments"))
\S4method{model.matrix}{snonlinearModel}(object,
type=c("regressors","instruments"))
}
\arguments{
\item{object}{Object of class \code{linearModel}, \code{rlinearModel} or any
system of equations class.}
\item{type}{Should the function returns the matrix of instruments or the
matrix of regressors. For \code{nonlinearModel} classes, type='regressors'
will produce an error message, because there is no such model matrix in
this case, at least not for now.}
}
\section{Methods}{
\describe{

\item{\code{signature(object = "linearModel")}}{
Linear models with not restrictions.
}

\item{\code{signature(object = "nonlinearModel")}}{
Nonlinear models with not restrictions.
}

\item{\code{signature(object = "rlinearModel")}}{
  linear models with restrictions.
}
  
 \item{\code{signature(object = "slinearModel")}}{
   System of linear equations with no restrictions. 
 }

 \item{\code{signature(object = "rslinearModel")}}{
   System of linear equations with restrictions. 
 }

 \item{\code{signature(object = "rsnonlinearModel")}}{
   System of nonlinear equations with restrictions. 
 }
 
 \item{\code{signature(object = "snonlinearModel")}}{
   System of nonlinear equations with no restrictions. 
 }
  }}
  \examples{
data(simData)

## Unrestricted model
model1 <- momentModel(y~x1+x2+x3, ~x2+x3+z1+z2, data=simData)
model.matrix(model1)[1:3,]

## Restrictions change the response
R <- c("x2=2","x3+x1=3")
rmodel1 <- restModel(model1, R)
rmodel1
model.matrix(rmodel1)[1:3,]
}
\keyword{methods}
