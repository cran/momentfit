\name{rhoFct}

\alias{rhoFct}
\alias{rhoEL}
\alias{rhoET}
\alias{rhoEEL}
\alias{rhoREEL}
\alias{rhoHD}
\alias{rhoETHD}
\alias{rhoETEL}
	
\title{GEL objective functions}

\description{
Functions that returns the GEL function \eqn{\rho(g(\theta,x)'\lambda)}
and its derivatives. 
}
\usage{
rhoET(gmat, lambda, derive = 0, k = 1)

rhoETEL(gmat, lambda, derive = 0, k = 1)

rhoEL(gmat, lambda, derive = 0, k = 1)

rhoEEL(gmat, lambda, derive = 0, k = 1)

rhoREEL(gmat, lambda, derive = 0, k = 1)

rhoHD(gmat, lambda, derive = 0, k = 1)

rhoETHD(gmat, lambda, derive = 0, k = 1)
}
\arguments{
\item{gmat}{The \eqn{n \times q} matrix of moments}

\item{lambda}{The \eqn{q \times 1} vector of Lagrange multipliers.} 

\item{derive}{An integer which indicates which derivative to return}

\item{k}{A numeric scaling factor that is required when \code{"gmat"} is
a matrix of time series which require smoothing. The value depends on
the kernel and is automatically set when the \code{"gelModels"} is
created.}

}


\references{
Anatolyev, S. (2005), GMM, GEL, Serial Correlation, and Asymptotic Bias. \emph{Econometrica}, \bold{73}, 983-1002.

Kitamura, Yuichi (1997), Empirical Likelihood Methods With Weakly Dependent Processes.
\emph{The Annals of Statistics}, \bold{25}, 2084-2102.

Kitamura, Y. and Otsu, T. and Evdokimov, K. (2013), Robustness,
Infinitesimal Neighborhoods and Moment Restrictions.
\emph{Econometrica}, \bold{81}, 1185-1201.

Newey, W.K. and Smith, R.J. (2004), Higher Order Properties of GMM and 
Generalized Empirical Likelihood Estimators. \emph{Econometrica}, \bold{72}, 219-255.

Smith, R.J. (2011), GEL Criteria for Moment Condition Models.
\emph{Econometric Theory}, \bold{27}(6), 1192--1235.
}

\value{
It returns the vector \eqn{\rho(gmat \lambda)} when \code{derive=0},
\eqn{\rho'(gmat \lambda)} when \code{derive=1} and \eqn{\rho''(gmat
  \lambda)} when \code{derive=2}.
 }


