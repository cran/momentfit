\name{kclassfit}
\alias{kclassfit}
\alias{getK}
\title{
K-Class Estimation Method
}
\description{
It estimates \code{linearModel} objects using the K-Class method. It
includes LIML, Fuller, TSLS and OLS.
}
\usage{
kclassfit(object, k, type = c("LIML", "Fuller", "BTSLS"), alpha = 1)

getK(object, alpha = 1, returnRes = FALSE)
}

\arguments{
  \item{object}{
    A model of class \code{linearModel}
  }
  \item{k}{
    The numeric value of k for the K-Class estimator. If missing, the value
    for LIML or Fuller is used. 
  }
  \item{type}{
    Which value of \code{k} should we use to fit the model? Only used if
    \code{k} is missing. 
  }
  \item{alpha}{A parameter for the Fuller method}

  \item{returnRes}{Should the function return the matrix of residuals
    from the first stage regression?}
}

\details{
  Let the model be \eqn{Y=X\beta+U}{Y=XB+U} and the matrix of instruments
  be \eqn{Z}{Z}. The K-Class estimator is
  \eqn{(X'(I-kM_z)X)^{-1}(X'(I-kM_z)Y)}. The function \code{getK} can be
  used to compute the value of \code{k} for both LIML and Fuller. When
  \code{type="BTSLS"}, the bias-adjusted TSLS of Nagar (1959) is
  computed. 
}

\value{
  The function returns an object of class \code{kclassfit}.
}  
\examples{

data(simData)
theta <- c(beta0=1,beta1=2)
model1 <- momentModel(y~x1, ~z1+z2, data=simData)
kclassfit(model1, type="LIML")
}
 
\keyword{LIML}
\keyword{Fuller}
