\name{hypothesisTest-methods}
\docType{methods}
\alias{hypothesisTest}
\alias{hypothesisTest-methods}
\alias{hypothesisTest,gmmfit,gmmfit-method}
\alias{hypothesisTest,missing,gmmfit-method}
\alias{hypothesisTest,gmmfit,missing-method}
\alias{hypothesisTest,sgmmfit,missing-method}
\alias{hypothesisTest,sgmmfit,sgmmfit-method}
\alias{hypothesisTest,missing,sgmmfit-method}
\title{ ~~ Methods for Function \code{hypothesisTest} in Package \pkg{momentfit} ~~}
\description{
Performs hypothesis tests on the coefficients estimated by any GMM fit method.
}
\usage{
\S4method{hypothesisTest}{gmmfit,missing}(object.u, object.r, R,
rhs=NULL, vcov=NULL, \dots)

\S4method{hypothesisTest}{sgmmfit,missing}(object.u, object.r, R,
rhs=NULL, vcov=NULL, \dots)

\S4method{hypothesisTest}{gmmfit,gmmfit}(object.u, object.r,
type=c("Wald", "LR", "LM"), sameVcov=TRUE, vcov=NULL,
firstStepWeight=FALSE, wObj=NULL, \dots)

\S4method{hypothesisTest}{sgmmfit,sgmmfit}(object.u, object.r,
type=c("Wald", "LR", "LM"), sameVcov=TRUE, vcov=NULL,
firstStepWeight=FALSE, wObj=NULL, \dots)

\S4method{hypothesisTest}{missing,gmmfit}(object.u, object.r, wObj=NULL)

\S4method{hypothesisTest}{missing,sgmmfit}(object.u, object.r, wObj=NULL)
}

\arguments{
  \item{object.u}{An object of class \code{gmmfit} or \code{sgmmfit}
    obtained using an unrestricted \code{"momentModel"} or
    \code{"sysModel"}.}
  \item{object.r}{An object of class \code{gmmfit} obtained using a
      restricted \code{"momentModel"} or \code{"sysModel"}.}
  \item{R}{If it is an object of class \code{gmmfit}, one of the model
    fit must be the restricted version of the other. The restrictions are
    then tested. If \code{R} is a \code{character} type, it expresses
    the restrictions using the coefficient names. If it \code{numeric}, it
    must be a matrix and the restrictions are \eqn{R\theta=0} for
    \code{NULL} rhs, or \eqn{R\theta=rhs} otherwise. If missing, the
    \code{gmmfit} must be a fitted restricted model, in which case, a LM
    test is performed.}
  \item{rhs}{A vector of right hand sides if \code{R} is \code{numeric}}
  \item{type}{Should we perform a Wald,  LR or LM test?}
  \item{sameVcov}{For the LR test, should we use the same estimate of
    the covariance matrix of the moment conditions? See details below.}
  \item{vcov}{For the Wald test, it is possible to provide the method
    with the covariance matrix of the coefficients.}
  \item{wObj}{For the LR test, it is possible to provide the
    \code{gmmWeights} object. In that case, the provided gmm weights
    object if used for the restricted and unrestricted models.}
  \item{\dots}{Other argument to pass to \code{\link{specTest}}.}
  \item{firstStepWeight}{Should we use the first step weighting matrix
  to compute the test (By default, the optimal weighting matrix is
  recomputed using the final vector of coefficient estimates). See
  details below.}
}
\section{Methods}{
\describe{

\item{\code{signature(object.u = "gmmfit", object.r = "gmmfit")}}{
Used to test a restricted model against an unrestricted one.
}

\item{\code{signature(object.u = "sgmmfit", object.r = "sgmmfit")}}{
Used to test a restricted model against an unrestricted one (for systems
of equations).
}

\item{\code{signature(object.u = "missing",  object.r= "gmmfit")}}{
Used to test a restricted model using the LM test.
}

\item{\code{signature(object.u = "missing",  object.r= "sgmmfit")}}{
Used to test a restricted model using the LM test (for systems of
equations).
}

\item{\code{signature(object.u = "gmmfit", object.r = "missing")}}{
Perform a Wald test using an unrestricted model and a restriction matrix
or vector.
}

\item{\code{signature(object.u = "sgmmfit", object.r = "missing")}}{
Perform a Wald test using an unrestricted model and a restriction matrix
or vector in systems of linear equations.
}
}}
\details{
  The LR test is the difference between the J-tests of the restricted
  and unrestricted models. It is therefore \eqn{n\bar{g}_r'W_r\bar{g}_r
    - n\bar{g}_u'W_u\bar{g}_u}, where \eqn{\bar{g}_r} and
  \eqn{\bar{g}_u} are respectively the restricted and unrestricted
  sample mean of the moment conditions, and \eqn{W_r} and \eqn{W_u}
  their respective optimal weigthing matrix. The test is therefore
  invalid if either of the weighting matrices does not converge to the
  inverse of the covariance matrix of the moment conditions. The
  restricted and unrestricted models must therefore be estimated by
  efficient GMM. This is not required for the Wald test.

  Asymptotically, it makes no difference which consistent estimate of
  \eqn{W_u} or \eqn{W_r} is used. However, it will make a difference in
  finite samples.

  If \code{sameVcov=TRUE}, both \eqn{W_r} and \eqn{W_u} are equal to the
  the optimal weighting matrix from the unrestricted model if
  \code{firstStepWeight=FALSE}, and they are equal to the first step
  weighting matrix (or the last step for iteratice GMM) of the
  unrestricted model if it is \code{TRUE}. For CUE, the value of
  \code{firstStepWeight} makes no difference since the weighting matrix
  and coefficients are computed simultaneously. Having \eqn{W_r=W_u}
  prevents the test to be negative in small samples.

  If \code{wObj} is provided, both \eqn{W_r} and \eqn{W_u} are equal to
  it. Of cource, \code{wObj} must be a consistent estimate of the
  optimal weighting matrix for the test to be valid.
}
\examples{
data(simData)

## Unrestricted model
model1 <- momentModel(y~x1+x2+x3, ~x2+x3+z1+z2+z3, data=simData, vcov="MDS")
res1 <- gmmFit(model1)

## Wald test 
R <- c("x1=0.5","x2=x3")
hypothesisTest(object.u=res1, R=R)

## LR tests

rmodel1 <- restModel(model1, R)
res2 <- gmmFit(rmodel1)
hypothesisTest(object.u=res1, object.r=res2, type="LR")

### LR and Wald should be the same as long as the same weighting
### matrix if used for both GMM fits, for the LR and Wald as well

# Unrestricted model and save the weights
res1 <- gmmFit(model1)
w <- res1@wObj
# estimate models with the same weights
res2 <- gmmFit(rmodel1, weights=w)

# LR test with the same weights
hypothesisTest(res1, res2, type="LR", wObj=w)

# Wald test with vcov based on the same weights (or the bread)
hypothesisTest(object.u=res1, R=R, breadOnly=TRUE)

### Another example with real data
data(Mroz)
model <- momentModel(log(wage)~educ+exper+I(exper^2),
                  ~exper+I(exper^2)+fatheduc+motheduc, vcov="MDS",
                  data=Mroz, centeredVcov=FALSE)
R <- c("educ=0","I(exper^2)=0")
rmodel <- restModel(model, R)

res1 <- gmmFit(model)
res2 <- gmmFit(rmodel, weights=res1@wObj)

hypothesisTest(object.u=res1, object.r=res2, type="LR", wObj=res1@wObj)
hypothesisTest(object.u=res1, object.r=res2, type="Wald",
vcov=vcov(res1, breadOnly=TRUE))

## LM test (identical to the other two tests as well)

hypothesisTest(object.r=res2)
# or 
hypothesisTest(object.u=res1, object.r=res2, type="LM")

## Wald with the Delta Method:
## To impose nonlinear restrictions, we need to convert
## the linear model into a nonlinear one
NLmodel <- as(model1, "nonlinearModel")
R1 <- c("theta2=2", "theta3=theta4^2")
rNLmodel <- restModel(NLmodel, R1)
res.u <- gmmFit(NLmodel)
res.r <- gmmFit(rNLmodel)
hypothesisTest(object.u=res.u, R=R1)

## LM

hypothesisTest(object.r=res.r)

## LR

hypothesisTest(object.r=res.r, object.u=res.u, type="LR")




}
\keyword{methods}
