\name{momentModel}

\alias{momentModel}
	
\title{Constructor for \code{"momentModel"} classes}

\description{
It builds an object class \code{"momentModel"}, which is a
union class for \code{"linearModel"},  \code{"nonlinearModel"},
\code{"formulaModel"} and \code{"functionModel"} classes. These are
classes for moment based models. This is the first step before running
any estimation algorithm.
}
\usage{
momentModel(g, x=NULL, theta0=NULL,grad=NULL,
            vcov = c("iid", "HAC", "MDS", "CL"),
            vcovOptions=list(), centeredVcov = TRUE, data=parent.frame(),
            na.action="na.omit", survOptions=list(), smooth=FALSE)
}
\arguments{
\item{g}{A function of the form \eqn{g(\theta,x)} and which returns a
  \eqn{n \times q} matrix with typical element \eqn{g_i(\theta,x_t)} for
  \eqn{i=1,...q} and \eqn{t=1,...,n}. This matrix is then used to build
  the q sample moment conditions. It can also be a formula if the model
  is linear (see detailsbelow).}

\item{x}{The matrix or vector of data from which the function
  \eqn{g(\theta,x)} is computed. If "g" is a formula, it is an \eqn{n
    \times Nh} matrix of instruments or a formula (see details below).}

\item{theta0}{A \eqn{k \times 1} vector of starting values. It is required
  only when "g" is a function because only then a numerical algorithm is
  used to minimize the objective function. If the dimension of
  \eqn{\theta} is one, see the argument "optfct".}

\item{grad}{A function of the form \eqn{G(\theta,x)} which returns a
  \eqn{q\times k} matrix of derivatives of \eqn{\bar{g}(\theta)} with
  respect to \eqn{\theta}. By default, the numerical algorithm
  \code{numericDeriv} is used. It is of course strongly suggested to
  provide this function when it is possible. This gradient is used to
  compute the asymptotic covariance matrix of \eqn{\hat{\theta}} and to
  obtain the analytical gradient of the objective function if the method
  is set to "CG" or "BFGS" in \code{\link{optim}} and if "type" is not
  set to "cue". If "g" is a formula, the gradiant is not required (see
  the details below).}

\item{vcov}{Assumption on the properties of the moment conditions. By
  default, they are weakly dependant processes. For \code{MDS}, we
  assume that the conditions are martingale difference sequences, which
  implies they are serially uncorrelated, but may be
  heteroscedastic. There is a difference between \code{iid} and
  \code{MDS} only when \code{g} is a formula. In that case, residuals
  are assumed homoscedastic as well as serially uncorrelated. For type
  \code{CL}, clustered covariance matrix is computed. The options are
  then included in \code{vcovOptions} (see \code{\link{meatCL}}). }

\item{vcovOptions}{A list of options for the covariance matrix of the
  moment conditions. See \code{\link{vcovHAC}} for the default values.}

\item{centeredVcov}{Should the moment function be centered when
  computing its covariance matrix. Doing so may improve inference.}

\item{data}{A data.frame or a matrix with column names (Optional). }

\item{na.action}{Action to take for missing values. If missing values
  are present and the option is set to \code{"na.pass"}, the model won't
  be estimable.}

\item{survOptions}{If needed, a list with the type of survey weights and
  the weights as a numeric vector, data.frame or formula. The type is either
  \code{"sampling"} or \code{"fequency"}.}

\item{smooth}{If \code{TRUE}, the moment function is smoothed using a
  kernel method.}
}

\value{
'momentModel' returns an object of one of the subclasses of \code{"momentModel"}.
 }

\references{
 Andrews DWK (1991),
  Heteroskedasticity and Autocorrelation Consistent Covariance Matrix Estimation.
  \emph{Econometrica}, \bold{59},
  817--858.

 Newey WK & West KD (1987), A Simple, Positive Semi-Definite,
 Heteroskedasticity and Autocorrelation Consistent Covariance
 Matrix. \emph{Econometrica}, \bold{55}, 703--708.

 Newey WK & West KD (1994), Automatic Lag Selection in Covariance
 Matrix Estimation. \emph{Review of Economic Studies}, \bold{61}, 631-653.
}
\examples{
data(simData)
theta <- c(beta0=1,beta1=2)

## A linearModel
model1 <- momentModel(y~x1, ~z1+z2, data=simData)

## A nonlinearModel
g <- y~beta0+x1^beta1
h <- ~z1+z2
model2 <- momentModel(g, h, c(beta0=1, beta1=2), data=simData)

## A functionModel
fct <- function(tet, x)
    {
        m1 <- (tet[1] - x)
        m2 <- (tet[2]^2 - (x - tet[1])^2)
        m3 <- x^3 - tet[1]*(tet[1]^2 + 3*tet[2]^2)
        f <- cbind(m1, m2, m3)
        return(f)
    }
dfct <- function(tet, x)
        {
        jacobian <- matrix(c( 1, 2*(-tet[1]+mean(x)), -3*tet[1]^2-3*tet[2]^2,0, 2*tet[2],
			   -6*tet[1]*tet[2]), nrow=3,ncol=2)
        return(jacobian)
        }
model3 <- momentModel(fct, simData$x3, theta0=c(beta0=1, beta1=2), grad=dfct)
}

