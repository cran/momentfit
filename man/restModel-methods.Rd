\name{restModel-methods}
\docType{methods}
\alias{restModel}
\alias{restModel-methods}
\alias{restModel,linearModel-method}
\alias{restModel,formulaModel-method}
\alias{restModel,slinearModel-method}
\alias{restModel,snonlinearModel-method}
\alias{restModel,nonlinearModel-method}
\alias{restModel,functionModel-method}

\title{ ~~ Methods for Function \code{restModel} in Package \pkg{momentfit} ~~}
\description{
  It converts \code{momentModel} objects into its restricted
  counterpart.
}
\usage{
\S4method{restModel}{linearModel}(object, R, rhs=NULL)

\S4method{restModel}{slinearModel}(object, R, rhs=NULL)

\S4method{restModel}{snonlinearModel}(object, R, rhs=NULL)

\S4method{restModel}{nonlinearModel}(object, R, rhs=NULL)

\S4method{restModel}{formulaModel}(object, R, rhs=NULL)

\S4method{restModel}{functionModel}(object, R, rhs=NULL)

}
\arguments{
  \item{object}{An object of class \code{"momentModel"} or
    \code{"sysModel"}.}  
  \item{R}{Either a matrix or a vector of characters for
    linear models and a list of formulas for nonlinear models}
  \item{rhs}{The right hand side of the linear restrictions. It is
    ignored for nonlinear models.}
}
\section{Methods}{
\describe{

\item{\code{signature(object = "linearModel")}}{
Method for object of class \code{linearModel}.
}

\item{\code{signature(object = "linearGel")}}{
  Method for all classes related to \code{linearGel}.
}

\item{\code{signature(object = "slinearModel")}}{
Method for object of class \code{slinearModel}.
}

\item{\code{signature(object = "snonlinearModel")}}{
Method for object of class \code{snonlinearModel}.
}

\item{\code{signature(object = "nonlinearModel")}}{
Method for object of class \code{nonlinearModel}.
}

\item{\code{signature(object = "nonlinearGel")}}{
Method for object of class \code{nonlinearGel}.
}

\item{\code{signature(object = "functionModel")}}{
Method for object of class \code{functionModel}.
}

\item{\code{signature(object = "functionGel")}}{
Method for object of class \code{functionGel}.
}

\item{\code{signature(object = "formulaModel")}}{
Method for object of class \code{formulaModel}.
}

\item{\code{signature(object = "formulaGel")}}{
Method for object of class \code{formulaGel}.
}
}}
\examples{
data(simData)
theta <- c(beta0=1,beta1=2)

## Unrestricted model
model1 <- momentModel(y~x1+x2+x3+z1, ~x1+x2+z1+z2+z3+z4, data=simData)

## Using matrix R
R <- matrix(c(1,1,0,0,0,0,0,2,0,0,0,0,0,1,-1),3,5, byrow=TRUE)
q <- c(0,1,3)

rmodel1 <- restModel(model1, R, q)
rmodel1

## Using character
## Many ways to write the constraints

R1 <- c("x1","2*x2+z1=2", "4+x3*5=3")
rmodel1 <- restModel(model1, R1)
rmodel1

## Works with interaction and identity function I()

model1 <- momentModel(y~x1*x2+exp(x3)+I(z1^2), ~x1+x2+z1+z2+z3+z4, data=simData)
R1 <- c("x1","exp(x3)+2*x1:x2", "I(z1^2)=3")
rmodel1 <- restModel(model1, R1)
rmodel1

## nonlinear constraints on a linear model
## we need to convert the linear model into a nonlinear one

model <- momentModel(y~x1+x2+x3+z1, ~x1+x2+z1+z2+z3+z4, data=simData)
NLmodel <- as(model, "nonlinearModel")

## To avoid having unconventional parameter names, which happens
## when I() is used or with interaction, the X's and coefficients are
## renamed

NLmodel@parNames

## Restriction can be a list of formula or vector of characters
## For the latter, it will be converted into a list of formulas

R1 <- c("theta2=2", "theta3=theta4^2")
rmod1 <- restModel(NLmodel, R1)
res1 <- gmmFit(rmod1)
res1
## recover the orignial form
coef(rmod1, coef(res1))

## with formulas

R2 <- list(theta2~2, theta3~1/theta4)
rmod2 <- restModel(NLmodel, R2)
res2 <- gmmFit(rmod2)
res2
coef(rmod2, coef(res2))

## The same can be done with function based models

}

\keyword{methods}
