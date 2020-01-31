\name{gel4}
\alias{gel4}

\title{GEL estimation}
\description{
The main functions and methods to fit any model with GEL. As opposed to
\code{\link{gelFit}}, models don't need to be created. It is all done by
the functions. It is meant to be more user friendly. 
}
\usage{

gel4(g, x=NULL, theta0=NULL,lambda0=NULL, getVcov=FALSE, 
     gelType = c("EL","ET","EEL","HD", "REEL","ETEL","ETHD"),
     vcov = c("MDS","iid","HAC"), grad=NULL,
     vcovOptions=list(), centeredVcov = TRUE,
     cstLHS=NULL, cstRHS=NULL, lamSlv=NULL,
     rhoFct=NULL, initTheta=c("gmm", "theta0"),
     data = parent.frame(),
     coefSlv=c("optim","nlminb","constrOptim"),
     smooth=FALSE, 
     lControl=list(), tControl=list())

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
  only when "g" is a function, a formula or a list of formulas. For
  these cases, they are needed to create the \code{"momentModel"} object.}

\item{lambda0}{The \eqn{q \times 1} vector of starting values for the
  Lagrange multipliers. By default a zero vector is used.}

\item{getVcov}{Should the method computes the covariance matrices of the
  coefficients and Lagrange multipliers.}

\item{gelType}{A character string specifying the type of GEL.}

\item{vcov}{Assumption on the properties of the moment conditions.}

\item{grad}{A function of the form \eqn{G(\theta,x)} which returns a
  \eqn{q\times k} matrix of derivatives of \eqn{\bar{g}(\theta)} with
  respect to \eqn{\theta}.}

\item{vcovOptions}{A list of options for the covariance matrix of the
  moment conditions. See \code{\link{vcovHAC}} for the default values.}

\item{centeredVcov}{Should the moment function be centered when
  computing its covariance matrix. Doing so may improve inference.}

\item{cstLHS}{The left hand side of the constraints to impose on the
  coefficients. See \code{\link{restModel}} for more details.}

\item{cstRHS}{The right hand side of the constraints to impose on the
  coefficients. See \code{\link{restModel}} for more details.}

\item{lamSlv}{An alternative solver for the Lagrange multiplier. By
  default, either \code{\link{Wu_lam}}, \code{\link{EEL_lam}},
  \code{\link{REEL_lam}} or \code{\link{getLambda}} is
  used. See the vignette for the required format.}

\item{rhoFct}{An optional function that return \eqn{\rho(v)}. This is
  for users who want a GEL model that is not built in the package. The
  four arguments of the function must be \code{"gmat"}, the matrix of
  moments, \code{"lambda"}, the vector of Lagrange multipliers,
  \code{"derive"}, which specify the order of derivative to return, and
  \code{k} a numeric scale factor required for time series and kernel
  smoothed moments.}

\item{initTheta}{Method to obtain the starting values for the
  coefficient vector. By default the GMM estimate with identity matrix
  is used. The second argument means that \code{"theta0"} is used
  instead. }

\item{data}{A required data.frame, in which all variables in g and x can
  be found.}

\item{smooth}{If TRUE, \code{"vcov"} is set to \code{"MDS"} and the
  moment conditions are smoothed using a kernel. See the vignette for
  more details.}

  \item{coefSlv}{Minimization solver for the coefficient vector.}

  \item{lControl}{A list of controls for the Lagrange multiplier
    algorithm.}

  \item{tControl}{A list of controls for the coefficient algorithm.}
}

\value{
It returns an object of class \code{"gelfit"}  
}
\references{
Anatolyev, S. (2005), GMM, GEL, Serial Correlation, and Asymptotic Bias. \emph{Econometrica}, \bold{73}, 983-1002.

Andrews DWK (1991), Heteroskedasticity and Autocorrelation Consistent Covariance Matrix Estimation.
\emph{Econometrica}, \bold{59}, 817--858.

Kitamura, Yuichi (1997), Empirical Likelihood Methods With Weakly Dependent Processes.
\emph{The Annals of Statistics}, \bold{25}, 2084-2102.

Kitamura, Y. and Otsu, T. and Evdokimov, K. (2013), Robustness,
Infinitesimal Neighborhoods and Moment Restrictions.
\emph{Econometrica}, \bold{81}, 1185-1201.

Newey, W.K. and Smith, R.J. (2004), Higher Order Properties of GMM and 
Generalized Empirical Likelihood Estimators. \emph{Econometrica}, \bold{72}, 219-255.

Smith, R.J. (2004), GEL Criteria for Moment Condition Models. \emph{Working paper, CEMMAP}.
}

\seealso{
\code{\link{gelFit}}, \code{\link{momentModel}} 
}
\examples{
data(simData)
res <- gel4(y~x1, ~z1+z2, vcov="MDS", gelType="ET", data=simData)
res
}
\keyword{GEL}
\keyword{Instruments}
\keyword{Moment}

