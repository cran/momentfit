\name{kernapply-methods}
\docType{methods}
\alias{kernapply}
\alias{kernapply-methods}
\alias{kernapply,momentModel-method}

\title{A kernel smoothing utility for \code{"momentModel"} classes}

\description{
It either generates the optimal bandwidth and kernel
weights or the smoothed moments of moment based models.
}

\usage{

\S4method{kernapply}{momentModel}(x, theta=NULL, smooth=TRUE, \dots)

}

\arguments{
  \item{x}{An object of class \code{"momentModel"}.}
  
  \item{theta}{An optional vector of coefficients. For
    \code{smooth=FALSE}, it is used to obtain the optimal bandwidth. If
    \code{NULL}, the bandwidth is obtained using one step GMM with the
    identity matrix as weights. For \code{smooth=TRUE}, the coefficient is
    required since the function returns the smoothed moments at a given
    vector of coefficients.}
  \item{smooth}{By default, it returns the smoothed moment matrix. If
    \code{FALSE}, it computes the optimal bandwidth and kernel weights.}
  \item{\dots}{Other arguments to pass. Currently not used}
}
  
\value{
  It return an object of class \code{"sSpec"}.
}

\references{ Anatolyev, S. (2005), GMM, GEL, Serial Correlation, and
Asymptotic Bias. \emph{Econometrica}, \bold{73}, 983-1002.

Kitamura, Yuichi (1997), Empirical Likelihood Methods With Weakly Dependent Processes.
\emph{The Annals of Statistics}, \bold{25}, 2084-2102.

Smith, R.J. (2011), GEL Criteria for Moment Condition Models.
\emph{Econometric Theory}, \bold{27}(6), 1192--1235.
}

\examples{
data(simData)
theta <- c(beta0=1,beta1=2)

## A linearModel
model1 <- momentModel(y~x1, ~z1+z2, data=simData,vcov="HAC",vcovOptions=list(kernel="Bartlett"))

### get the bandwidth
### Notice that the kernel name is the not the same
### That's because a Truncated kernel for smoothing
### lead to a Bartlett kernel for the HAC of the moments
### See Smith (2011)
kernapply(model1, smooth=FALSE)


### Adding the kernel option to the model

model2 <- momentModel(y~x1, ~z1+z2,
data=simData,vcov="HAC",vcovOptions=list(kernel="Bartlett"), smooth=TRUE)

kernapply(model2, theta)$smoothx[1:5,]

}

