\name{weakTest}
\alias{CDtest}
\alias{SWtest}
\alias{MOPtest}
\alias{SYTables}
\alias{LewMertest}
\title{
Tests for weak Instruments
}
\description{
This is a collection of tests for weak instruments. It includes the
Cragg and Donald test for the reduced rank hypothesis, the conditional
F-test of Sanderson and Windmeijer (2016), the robust effective F
test of Montiel Olea and Pflueger (2013) and the Lewis and Mertens
(2022) robust test with multiple endogenous variables. The critical
values of Stock and Yogo (2005) for the null hypothesis of weak
instruments are also provided.}

\usage{

SWtest(object, j=1, print=TRUE, ...)

CDtest(object, print=TRUE, SWcrit=FALSE, ...)

MOPtest(object, tau=0.10, size=0.05, print=TRUE,
        estMethod = c("TSLS", "LIML"), simplified = TRUE,
        digits = max(3L, getOption("digits") - 3L), ...)

LewMertest(object, tau=0.10, size=0.05, print=TRUE,
           simplified = TRUE, digits = max(3L, getOption("digits") -
           3L), npoints=10, ...)

SYTables(object, print=TRUE, SWcrit=FALSE)

}

\arguments{
  \item{object}{
    A model of class \code{linearModel}
  }

  \item{j}{The Sanderson and Windmeijer test is based on the regression
  of \eqn{y_j-\delta y_{-j}} on the set on instruments, where \eqn{y_j}
  is the jth included endogenous variable and \eqn{y_{-j}} is the vector
  of the other included endogenous variables.
  }

  \item{tau}{The desired bias.}

  \item{size}{The size of the test for the critical value.}

  \item{estMethod}{Which method we want to test the weak instruments
  for? The method determined the critical value associated with a given
  relative bias.}

\item{simplified}{Should we perform the simplified test? The test
  produce more conservative critical values. The generalized test is
  computationally more demanding.}
  
  \item{print}{
    If \code{TRUE}, the result is printed and nothing is returned. See
    below for details on what the functions return when it is set to
    \code{FALSE}
  }
  
  \item{digits}{The number of digits to print when \code{print} is set
  to \code{TRUE}.}

  \item{SWcrit}{If true, the critical values and the test are ajusted to
  Sanderson and Windmeijer (2016). See details.}

\item{npoints}{The number of minimizations to run in order to find the
  global minimum, which is used for the Lewis and Mertens critical
  value. It affects the results only when \code{simplified} is set to
  \code{TRUE}.} 

  \item{...}{Arguments passed to \code{\link{formatC}} to print the
  statistic. For \code{MOPtest}, they are passed to
  \code{momentfit:::getMOPx}. 
  }
}

\details{ The \code{CDtest} function computes the test for the model
  object and runs \code{SYTables} to extract the critical values that are
  specific to the model.

  Let \eqn{l} be the number of included endogenous variables and \eqn{s}
  be the number of excluded exogenous. The Stock and Yogo (2005) tables
  are based on these two values and they are only available for
  \eqn{l=1,2} and \eqn{s=1,...,30}. For the SW test, which is only
  defined for $l>1$, we compare the statistic with the critical values
  associated with \eqn{l=l-1} and \eqn{s=s-l+1}. These are the critical
  values that are returned when the argument \code{SWcrit} of
  \code{SYTables} is set to \code{TRUE}. This allows us to test for weak
  instruments in models with 2 or 3 included endogenous variable using
  the same tables of critical values.

  Sanderson and Windmeijer (2016) show that we can use the same critical
  values if we modify the CD test by multiplying it by \eqn{s/(s-l+1)},
  which is what the function \code{CDtest} returns if the argument
  \code{SWcrit} is set to \code{TRUE}.

  }

\value{
  The function \code{CDtest} returns the Cragg anb Donald statistic when
  \code{print} is set to \code{FALSE}.

  The function \code{SYTable} returns a list with the following elements
  when \code{print} is set to \code{FALSE}:

  \item{biasTSLS, biasFuller}{Named vectors of critical values for TSLS
  or Fuller (see \code{\link{kclassfit}}). These critical values are
  used to achieve the maximum relative bias given by their names. The
  values are defined only for models with a number of overidentifying
  restrictions greater or equal to 2.} 

  \item{sizeTSLS, sizeLIML}{A named vector of critical values for TSLS
  or LIML (see \code{\link{kclassfit}}). These critical values are used
  to acheive a size that does not exceed the one given by their names.}

 Both functions return \code{NULL} invisibly when \code{print} is set to
 \code{TRUE}

}  
\examples{

data(simData)
theta <- c(beta0=1,beta1=2)
model1 <- momentModel(y~x1, ~z1+z2, data=simData)
CDtest(model1)
model2 <- momentModel(y~y1+y2+x1, ~z1+z2+z3+z4+x1, data=simData)
SWtest(model2, 1, FALSE)
SWtest(model2, 2, FALSE)
}
 
\keyword{LIML}
\keyword{Fuller}
\keyword{Weak Instruments Test}
\keyword{Two Stage Least Squares}

