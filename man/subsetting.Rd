\name{[-methods}
\docType{methods}
\alias{[-methods}
\alias{subset}
\alias{subset,regModel-method}
\alias{subset,functionModel-method}
\alias{subset,formulaModel-method}
\alias{subset,sysModel-method}
\alias{[,momentWeights,numeric,missing-method}
\alias{[,momentWeights,missing,missing-method}
\alias{[,sysMomentWeights,numeric,list-method}
\alias{[,sysMomentWeights,numeric,missing-method}
\alias{[,sysMomentWeights,missing,list-method}
\alias{[,momentModel,missing,missing-method}
\alias{[,regModel,numeric,missing-method}
\alias{[,functionModel,numeric,missing-method}
\alias{[,formulaModel,numeric,missing-method}
\alias{[,rfunctionModel,numeric,missing-method}
\alias{[,sysModel,missing,list-method}
\alias{[,snonlinearModel,numeric,missing-method}
\alias{[,sfunctionModel,numeric,missing-method}
\alias{[,slinearModel,numeric,missing-method}
\alias{[,rsnonlinearModel,numeric,missing-method}
\alias{[,rslinearModel,numeric,missing-method}
\alias{[,sysModel,numeric,list-method}
\alias{[,sysModel,missing,missing-method}

\title{Subsetting methods}
\description{
Different subsetting methods for S4 class objects of the package. The
subset method returns an new object with observations selected by the
second argument. See example.
}
\section{Methods}{
\describe{

\item{\code{signature(x = "momentWeights", i = "integer", j = "missing")}}{
  It creates a partition from the weighting matrix. 
}

\item{\code{signature(x = "momentWeights", i = "missing", j = "missing")}}{  
  It generates the whole weighting matrix. 
}


\item{\code{signature(x = "sysMomentWeights", i = "missing", j = "list")}}{
  It creates a partition from the weighting matrix. j has no effect
  here.  It creates a partition from the weighting matrix in a systemof  
  equations. i selects the equation and the list j the moments in each
  equation. Missing i means all equations.
}

\item{\code{signature(x = "sysMomentWeights", i = "numeric", j =
    "missing")}}{
    It creates a partition from the weighting matrix. j has no effect
  here.  It creates a partition from the weighting matrix in a systemof  
  equations. i selects the equation and the list j the moments in each
  equation. Missing j means all moments.
}

\item{\code{signature(x = "sysMomentWeights", i = "missing", j = "missing")}}{
  No effect. It returns x.
}


\item{\code{signature(x = "snonlinearModel", i = "numeric", j="missing")}}{
  It generates a system of equations with a subset of equations selected
  by i. If the number of remaining equations is one, it returns an
  object of class \code{"nonlinearGmm"}.
}


\item{\code{signature(x = "slinearModel", i = "numeric", j="missing")}}{
  It generates a system of equations with a subset of equations selected
  by i. If the number of remaining equations is one, it returns an
  object of class \code{"linearModel"}.
}

\item{\code{signature(x = "rslinearModel", i = "numeric", j="missing")}}{
  It is only use to select one equation when no cross-equation
  restrictions are imposed. Only one equation can be selected. 
}

\item{\code{signature(x = "rsnonlinearModel", i = "numeric", j="missing")}}{
  It is only use to select one equation when no cross-equation
  restrictions are imposed. Only one equation can be selected. 
}

\item{\code{signature(x = "sysMomentModel", i = "numeric", j="list")}}{
  It generates a system of equations with a subset of equations selected
  by i and a subset of moment conditions selected by j. If the number of
  remaining equations is one, it returns an object of class \code{"linearGmm"}.
}

\item{\code{signature(x = "sysMomentModel", i = "missing", j="missing")}}{
  No effect. It returns x.
}


\item{\code{signature(x = "momentModel", i = "missing", j = "missing")}}{
Returns the model without any change.
}

\item{\code{signature(x = "functionModel", i = "numeric", j = "missing")}}{
  It generates the same model with a subset of the moment conditions.
}

\item{\code{signature(x = "formulaModel", i = "numeric", j = "missing")}}{
  It generates the same model with a subset of the moment conditions.
}

\item{\code{signature(x = "rfuncionModel", i = "numeric", j = "missing")}}{
  It generates the same model with a subset of the moment conditions. j
  has no effect here.
}

}}
\examples{
data(simData)
model1 <- momentModel(y~x1+x2, ~x2+x3+z1+z2+z3, data=simData, vcov="MDS")
w <- evalWeights(model1, theta=1:3)
w[]
w[1:3]

## A model with a subset of the instruments
model1[1:4]

## Selecting the observations:

subset(model1, simData[["x1"]]<3)
subset(model1, 1:25)

}

\keyword{methods}
