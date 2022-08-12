\name{gmm4}
\alias{gmm4}
\alias{tsls,formula-method}
\alias{tsls,list-method}
\alias{ThreeSLS,list-method}

\title{GMM estimation
}
\description{
The main functions and methods to fit any model with GMM. As opposed to
\code{\link{gmmFit}}, models don't need to be created. It is all done by
the functions. It is meant to be more user friendly. This document
needs to changed. It is just a copy and paste from the gmm package
}
\usage{
gmm4(g, x, theta0 = NULL, grad = NULL, 
     type = c("twostep", "iter", "cue", "onestep"),
     vcov = c("iid", "HAC", "MDS", "TrueFixed", "CL"),
     initW = c("ident", "tsls", "EbyE"), weights = "optimal", 
     itermaxit = 50, cstLHS=NULL, cstRHS=NULL,
     vcovOptions=list(), survOptions=list(),
     itertol = 1e-07, centeredVcov = TRUE,
     data = parent.frame(), ...)

\S4method{tsls}{formula}(model, x, vcov = c("iid", "HAC", "MDS", "CL"),
         vcovOptions=list(), survOptions=list(), centeredVcov = TRUE,
         data = parent.frame())

\S4method{tsls}{list}(model, x=NULL, vcov = c("iid", "HAC", "MDS",
          "CL"), vcovOptions=list(), survOptions=list(),
          centeredVcov = TRUE, data = parent.frame())

\S4method{ThreeSLS}{list}(model, x=NULL, vcov = c("iid", "HAC", "MDS",
          "CL"), vcovOptions=list(), survOptions=list(),
          centeredVcov = TRUE, data = parent.frame())
}
\arguments{
\item{model}{A formula or a list of formulas.}  
\item{g}{A function of the form \eqn{g(\theta,x)} and which returns a
  \eqn{n \times q} matrix with typical element \eqn{g_i(\theta,x_t)} for
  \eqn{i=1,...q} and \eqn{t=1,...,n}. This matrix is then used to build
  the q sample moment conditions. It can also be a formula if the model
  is linear or nonlinear, or a list of formulas for systems of equations.}

\item{x}{The matrix or vector of data from which the function
  \eqn{g(\theta,x)} is computed. If "g" is a formula, it is an \eqn{n
    \times Nh} matrix of instruments or a formula (see details below).}

\item{theta0}{A \eqn{k \times 1} vector of starting values. It is required
  only when "g" is a function or a nonlinear equation defined by a
  formula, in which case, it must be a named vector}

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

  \item{type}{What GMM methods should we use? for
    \code{type=="onestep"}, if \code{"weights"} is not a matrix, the
    model will be estimated with the weights equals to the identity
    matrix}

  \item{vcov}{Assumption on the properties of the random vector x. By
  default, x is a weakly dependant process. The "iid" option will avoid
  using the HAC matrix which will accelerate the estimation if one is
  ready to make that assumption. The option "TrueFixed" is used only
  when the matrix of weights is provided and it is the optimal one. For
  type \code{CL}, clustered covariance matrix is computed. The options
  are then included in \code{vcovOptions} (see \code{\link{meatCL}}).}

  \item{initW}{How should be compute the initial coefficient vector in
    the first. It only makes a difference for linear models for which
    the choice is GMM with identity matrix or two-stage least quares.}

  \item{weights}{What weighting matrix to use? The
    choices are \code{"optimal"}, in which case it is the inverse of the
    moment vovariance matrix, \code{"ident"} for the identity matrix, or
    a fixed matrix.}
  
  \item{itermaxit}{Maximum iterations for iterative GMM}
  
  \item{itertol}{Tolance for the stopping rule in iterative GMM}

  \item{centeredVcov}{Should the moment function be centered when
    computing its covariance matrix. Doing so may improve inference.}
  
  \item{data}{A data.frame or a matrix with column names (Optional). }
  
  \item{cstLHS}{The left hand side of the constraints to impose on the
    coefficients. See \code{\link{restModel}} for more details.}
  
  \item{cstRHS}{The right hand side of the constraints to impose on the
    coefficients. See \code{\link{restModel}} for more details.}

  \item{vcovOptions}{A list of options for the covariance matrix of the
    moment conditions. See \code{\link{vcovHAC}} for the default values.}

  \item{survOptions}{If needed, a list with the type of survey weights and
    the weights as a numeric vector, data.frame or formula. The type is either
    \code{"sampling"} or \code{"fequency"}.}  
 
  \item{\dots}{Arguments to pass to \code{\link{optim}} when the model
    is nonlinear.}
}

\value{
It returns an object of class \code{"gmmfit"}  
}
\references{
  Zeileis A (2006), Object-oriented Computation of Sandwich Estimators.
  \emph{Journal of Statistical Software}, \bold{16}(9), 1--16.
  URL \doi{10.18637/jss.v016.i09}.

 Andrews DWK (1991),
  Heteroskedasticity and Autocorrelation Consistent Covariance Matrix Estimation.
  \emph{Econometrica}, \bold{59},
  817--858.

 Newey WK & West KD (1987), A Simple, Positive Semi-Definite,
 Heteroskedasticity and Autocorrelation Consistent Covariance
 Matrix. \emph{Econometrica}, \bold{55}, 703--708.

 Newey WK & West KD (1994), Automatic Lag Selection in Covariance
 Matrix Estimation. \emph{Review of Economic Studies}, \bold{61}, 631-653.

  Hansen, L.P. (1982),
  Large Sample Properties of Generalized Method of Moments Estimators.
  \emph{Econometrica}, \bold{50},
  1029-1054,

  Hansen, L.P. and Heaton, J. and Yaron, A.(1996),
  Finite-Sample Properties of Some Alternative GMM Estimators.
  \emph{Journal of Business and Economic Statistics}, \bold{14}
  262-280.
}

\seealso{
\code{\link{gmmFit}}, \code{\link{momentModel}} 
}
\examples{
data(simData)

res <- gmm4(y~x1, ~z1+z2, vcov="MDS", type="iter", data=simData)
res
}
\keyword{GMM}% use one of  RShowDoc("KEYWORDS")
\keyword{Instruments}% use one of  RShowDoc("KEYWORDS"
\keyword{Moment}% use one of  RShowDoc("KEYWORDS"

