\name{lambdaAlgo}

\alias{Wu_lam}
\alias{EEL_lam}
\alias{REEL_lam}
\alias{ETXX_lam}
\alias{getLambda}
	
\title{Algorithms to solve for the Lagrange multiplier}

\description{
The algorithms finds the vector or Lagrange multipliers that maximizes
the GEL objective function for a given vector of coefficient \eqn{\theta}.
}
\usage{
Wu_lam(gmat, tol=1e-8, maxiter=50, k=1)

EEL_lam(gmat, k=1) 

REEL_lam(gmat, tol=NULL, maxiter=50, k=1)

ETXX_lam(gmat, lambda0, k, gelType, algo, method, control)

getLambda(gmat, lambda0=NULL, gelType=NULL, rhoFct=NULL, 
          tol = 1e-07, maxiter = 100, k = 1, method="BFGS", 
          algo = c("nlminb", "optim", "Wu"), control = list(),
          restrictedLam=integer()) 
}
\arguments{
\item{gmat}{The \eqn{n \times q} matrix of moments}

\item{lambda0}{The \eqn{q \times 1} vector of starting values for the
  Lagrange multipliers.} 

\item{tol}{A tolerance level for the stopping rule in the Wu algorithm}

\item{maxiter}{The maximum number of iteration in the Wu algorithm}

\item{gelType}{A character string specifying the type of GEL. The
  available types are \code{"EL"}, \code{"ET"}, \code{"EEL"},
  \code{"HD"} and "REEL". For the latter, the algorithm restricts the
  implied probabilities to be non negative.} 

\item{rhoFct}{An optional function that return \eqn{\rho(v)}. This is
  for users who want a GEL model that is not built in the package. The
  four arguments of the function must be \code{"gmat"}, the matrix of
  moments, \code{"lambda"}, the vector of Lagrange multipliers,
  \code{"derive"}, which specify the order of derivative to return, and
  \code{k} a numeric scale factor required for time series and kernel
  smoothed moments.}

\item{k}{A numeric scaling factor that is required when \code{"gmat"} is
a matrix of time series which require smoothing. The value depends on
the kernel and is automatically set when the \code{"gelModels"} is
created.}

\item{method}{This is the method for \code{\link{optim}}.}

\item{algo}{Which algorithm should be used to maximize the GEL objective
  function. If set to \code{"Wu"}, which is only for \code{"EL"}, the Wu
  (2005) algorithm is used.}

\item{control}{A list of control to pass to either \code{\link{optim}}
  or \code{\link{nlminb}}.}

\item{restrictedLam}{A vector of integers indicating which
  \code{"lambda"} are restricted to be equal to 0.}
  
}

\details{The \code{ETXX_lam} is used for ETEL and ETHD. In general, it
 computes lambda using ET, and returns the value of the objective
 function determined by the \code{gelType}.  }

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

Wu, C. (2005), Algorithms and R codes for the pseudo empirical
likelihood method in survey sampling.
\emph{Survey Methodology}, \bold{31}(2), page 239.

}


\value{
It returns the vector \eqn{\rho(gmat \lambda)} when \code{derive=0},
\eqn{\rho'(gmat \lambda)} when \code{derive=1} and \eqn{\rho''(gmat
  \lambda)} when \code{derive=2}.
 }


