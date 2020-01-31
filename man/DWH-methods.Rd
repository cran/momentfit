\name{DWH-methods}
\docType{methods}
\alias{DWH}
\alias{DWH-methods}
\alias{DWH,gmmfit,lm-method}
\alias{DWH,gmmfit,missing-method}
\alias{DWH,gmmfit,gmmfit-method}
\title{ ~~ Methods for Function \code{DWH} in Package \pkg{momentfit} ~~}
\description{
It performs the Durbin-Wu-Hausman test on GMM fit models.
}
\usage{
\S4method{DWH}{gmmfit,missing}(object1, object2)

\S4method{DWH}{gmmfit,lm}(object1, object2,
tol=sqrt(.Machine$double.eps), v1=NULL, v2=NULL, \dots)

\S4method{DWH}{gmmfit,gmmfit}(object1, object2,
tol=sqrt(.Machine$double.eps), v1=NULL, v2=NULL, \dots)

}

\arguments{
  \item{object1}{Object of class \code{"gmmfit"}.}

  \item{object2}{Object of class \code{"gmmfit"} or \code{"lm"}. If
  missing, the DWH test is a two step test in which the fitted
  endogenous variables from the first step are added to the
  regression. In that case, the test a a test of significance of the
  coefficients of the fitted endogenous variables.}

  \item{v1}{Alternatively, we can provide a different covariance
    matrix for object1}

  \item{v2}{Alternatively, we can provide a different covariance
    matrix for object2}
    
  \item{tol}{Tolerance for the Moore-Penrose generalized inverse}
  
  \item{\dots}{Argument to pass to \code{\link{vcov}}}
}

\section{Methods}{
\describe{

\item{\code{signature(object1 = "gmmfit", object2 = "lm")}}{
%%  ~~describe this method here~~
}

\item{\code{signature(object1 = "gmmfit", object2 = "gmmfit")}}{
%%  ~~describe this method here~~
}

\item{\code{signature(object1 = "gmmfit", object2 = "missing")}}{
%%  ~~describe this method here~~
}
}}

\references{
Green, W.H.. (2012). \emph{Econometric Analysis, 7th edition}, Prentice Hall.    
}

\examples{

### Exampe 8.7 of Greene (2012)
data(ConsumptionG)
Y <- ConsumptionG$REALDPI
C <- ConsumptionG$REALCONS
n <- nrow(ConsumptionG)
Y1 <- Y[-n]; Y <- Y[-1]
C1 <- C[-n]; C <- C[-1]
dat <- data.frame(Y=Y,Y1=Y1,C=C,C1=C1)



model1 <- momentModel(C~Y, ~Y, data=dat, vcov="iid")
model2 <- momentModel(C~Y, ~Y1+C1, data=dat, vcov="iid")
res1 <- tsls(model1)
res2 <- tsls(model2)
res <- lm(C~Y)

## Exampke 8.7-2. The difference is explained by the rounding
## error in Greene. Only the first the 3 digits of the t-test are used.
DWH(res2)

## Example 8.7-1. Not quite the same.
DWH(res2, res1)

## using lm object to compare OLS and 2SLS:
## The same adjustment on the vcov must be done (it is by default in lm)
## otherwise the different in the covariance matrices is mostly caused by the
## different ways to compute them.
DWH(res2, res, df.adj=TRUE)

## To reproduce the same results as Exampke 8.7-1,
## we need to specify the variance.
## But it is not necessary as the above way is
## asymptotically equivalent
X <- model.matrix(model1)
Xhat <- qr.fitted(res2@wObj@w, X)
s2 <- sum(residuals(res)^2)/(res$df.residual)
v1 <-  solve(crossprod(Xhat))*s2
v2 <- solve(crossprod(X))*s2
DWH(res2, res, v1=v1, v2=v2)


}
\keyword{methods}

