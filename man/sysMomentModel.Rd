\name{sysMomentModel}

\alias{sysMomentModel}
	
\title{Constructor for \code{"sysMomentModel"} classes}

\description{
It builds the object of either class \code{"slinearModel"} or
\code{"snonlinearModel"}, which are system of equations based on moment
conditions.
}
\usage{
sysMomentModel(g, h=NULL, theta0=NULL, grad=NULL, 
            vcov = c("iid", "HAC", "MDS", "CL"),
            vcovOptions=list(), centeredVcov = TRUE,
            data=parent.frame(),na.action="na.omit",
            survOptions=list())
}
\arguments{
\item{g}{A list of linear or nonlinear regression formulas for each
  equation in the system.}

\item{h}{A list of linear formulas for the instruments in each
  equation in the system.}

\item{theta0}{A list of vectors of starting values. It is required
  only when the equations are nonlinear, in which case, it must be a
  list of named vector, with the names corresponding to the coefficient
  names in the regression formulas.}

\item{grad}{A list of functions that returns the derivative of the
  moment functions. Only used if \code{g} is a list of functions.}

\item{vcov}{Assumption on the properties of the moment conditions. By
  default, they are weakly dependant processes. For \code{MDS}, we
  assume that the conditions are martingale difference sequences, which
  implies they are serially uncorrelated, but may be
  heteroscedastic. There is a difference between \code{iid} and
  \code{MDS} only when \code{g} is a formula. In that case, residuals
  are assumed homoscedastic as well as serially uncorrelated. For type
  \code{CL}, clustered covariance matrix is computed. The options are
  then included in \code{vcovOptions} (see \code{\link{meatCL}}).}

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
}

\value{
'sysMomentModel' returns an object of one of the subclasses of
\code{"sysMomentModel"}.
 }

 \references{
 Hayashi, F. (2000). \emph{Econometrics}, New Jersey: Princeton
 University Press.
   
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
set.seed(1122)
x1 <- rchisq(50,5)
x2 <- rchisq(50,5)
x3 <- rnorm(50)
x4 <- rnorm(50)
z1 <- .2*x1+rnorm(50)
z2 <- .2*x2+rnorm(50)
z3 <- rnorm(50)
z4 <- rnorm(50)
z5 <- rnorm(50)
y1 <- x1+rnorm(50)
y2 <- 2*x1+rnorm(50)
y3 <- 0.5*x2+rnorm(50)
dat <- data.frame(y1=y1,y3=y3,y2=y2, z1=z1,x1=x1,z2=z2,x2=x2,z3=z3,x3=x3,
                  x4=x4,z4=z4,z5=z5)

g1 <- y1~x1+x4; h1 <- ~z1+z2+z3+z4+x4
g2 <- y2~x1+x2+x3; h2 <- ~z1+z2+z3+z4+x3
g3 <- y3~x2+x3+x4; h3 <- ~z2+z3+z4+x3+x4
g <- list(g1,g2,g3)
h <- list(h1,h2,h3)

smodel <- sysMomentModel(g, h, data=dat)

## not really nonlinear
nlg <- list(y1~theta0+theta1*x1+theta2*x4,
            y2~alpha0+alpha1*x1+alpha2*x2+alpha3*x3,
            y3~beta0+beta1*x2+beta2*x3+beta3*x4)
theta0 <- list(c(theta0=1,theta1=2,theta2=3),
              c(alpha0=1,alpha1=2,alpha2=3, alpha3=4),
              c(beta0=1,beta1=2,beta2=3,beta3=4))
snmodel <- sysMomentModel(nlg, h, theta0, data=dat)


}

