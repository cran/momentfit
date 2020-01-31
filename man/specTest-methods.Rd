\name{specTest-methods}
\docType{methods}
\alias{specTest}
\alias{specTest-methods}
\alias{specTest,gmmfit,missing-method}
\alias{specTest,gelfit,missing-method}
\alias{specTest,sgmmfit,missing-method}
\alias{specTest,gmmfit,numeric-method}
\title{ ~~ Methods for Function \code{specTest} in Package \pkg{momentfit} ~~}
\description{
It computes tests of specification for GMM fit.
}
\usage{
\S4method{specTest}{gmmfit,missing}(object, which, df.adj=FALSE, wObj=NULL)

\S4method{specTest}{sgmmfit,missing}(object, which, df.adj=FALSE, wObj=NULL)

\S4method{specTest}{gmmfit,numeric}(object, which)

\S4method{specTest}{gelfit,missing}(object, which,
                    type = c("All", "LR", "LM", "J"))
}
\arguments{
  \item{object}{GMM or GEL fit object}
  \item{which}{Which sub-moment conditions to test.}
  \item{df.adj}{Should we adjust the covariance matrix of the moment
    conditions for degrees of freedom. If \code{TRUE}
    the covariance matrix is multiplied by \code{n/(n-k)}, where
    \code{n} is the sample size and \code{k} is the number of
    coefficients. For heteroscedastic robust covariance matrix,
    adjusting is equivalent to computing HC1 while not adjusting is
    HC0.}
  \item{wObj}{An object of class \code{gmmWeights}. If NULL (the
    recommended value), the optimal weights is computed at the
    fitted coefficient estimates. It is used by
    \code{\link{hypothesisTest}} if one wants the LR statistics to be
    computed using the same weights for the restricted and unrestricted
    model.}
  \item{type}{For GEL, three specification tests are available}
}

\section{Methods}{
\describe{
\item{\code{signature(object = "gmmfit", which="missing")}}{
}

\item{\code{signature(object = "sgmmfit", which="missing")}}{
}

\item{\code{signature(object = "gmmfit", which="numeric")}}{
}
}}

\references{
Eichenbaum, M. and Hansen L. and Singleton, K. (1985). A time Series
Analysis of Representative Agent Models of Consumption and Leisure
Choise under Uncertainty. \emph{Quarterly Journal of Economics},
\bold{103}, 51--78.

Hayashi, F. (2000). \emph{Econometrics}, New Jersey: Princeton
University Press.
}

\examples{
data(simData)
model1 <- momentModel(y~x1, ~z1+z2, data=simData)

res <- gmmFit(model1)
specTest(res)

## Hayashi Example 3.3 (there is not result in the book but
## that's how we would do it for YEAR=1967
data(Griliches)
dat <- subset(Griliches, YEAR==67)
model <- momentModel(LW~S+EXPR+IQ, ~S+EXPR+AGE+MED, data=dat, vcov="MDS")
res <- gmmFit(model)
## testing the orthogonality conditions of S
specTest(res, 2)

}
\keyword{methods}
