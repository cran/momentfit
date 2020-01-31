\name{confint-methods}
\docType{methods}
\alias{confint-methods}
\alias{confint,ANY-method}
\alias{confint,gelfit-method}
\alias{confint,gmmfit-method}
\alias{confint,numeric-method}
\alias{confint,data.frame-method}
\alias{confint,matrix-method}
\alias{confint,ANY-method}
\title{ ~~ Methods for Function \code{confint} in Package \pkg{stats} ~~}
\description{
Method to contruct confidence intervals for objects of class
\code{"gmmfit"} and \code{"gelfit"}.
}

\usage{
\S4method{confint}{gmmfit}(object, parm, level = 0.95, vcov=NULL,
                    area=FALSE, npoints=50, \dots)

\S4method{confint}{gelfit}(object, parm, level = 0.95, lambda = FALSE,
                    type = c("Wald", "invLR", "invLM", "invJ"),
                    fact = 3, corr = NULL, vcov=NULL,
                    area = FALSE, npoints = 20, cores=4, \dots)

\S4method{confint}{numeric}(object, parm, level = 0.95, gelType="EL", 
                    type = c("Wald", "invLR", "invLM", "invJ"),
                    fact = 3, vcov="iid") 

\S4method{confint}{data.frame}(object, parm, level = 0.95, gelType="EL", 
                    type = c("Wald", "invLR", "invLM", "invJ"),
                    fact = 3, vcov="iid", npoints=10, 
                    cores=4) 

\S4method{confint}{matrix}(object, parm, level = 0.95, gelType="EL", 
                    type = c("Wald", "invLR", "invLM", "invJ"),
                    fact = 3, vcov="iid", npoints=10, 
                    cores=4) 

\S4method{confint}{ANY}(object, parm, level = 0.95, \dots)
}

\arguments{

  \item{object}{Object of class \code{"gmmfit"},
    \code{"gelfit"}, \code{"numeric"} or  \code{"data.frame"}.}

  \item{parm}{Vector of integers or characters for selecting the
    elements for which the intervals should be computed.}

  \item{level}{The confidence level.}

  \item{lambda}{Should be compute intervals for the Lagrange
    multipliers?}
  
  \item{type}{The type of confidence intervals. The default is the Wald
    interval, and the others are computed by inverting the LR, LM or J
    specification test.}

  \item{fact}{For the inversion of the specification tests,
    \code{\link{uniroot}} searches within \code{fact} standard error of
    the coefficient estimates}

  \item{corr}{Correction to apply to the specification tests}
  
  \item{vcov}{For Wald intervals, an optional covariance matrix can be
    provided. For \code{"numeric"} or \code{"data.frame"}, it specifies
    the type of observations.}

  \item{cores}{The number of cores for \code{\link{mclapply}}. It is set
    to 1 for Windows OS.}

  \item{gelType}{Type of GEL confidence interval.}

  \item{npoints}{Number of equally spaced points for the confidence
    region}

  \item{area}{If TRUE, a cnnfidence region is computed. The length of
  \code{"parm"} must be 2 in that case.}

  \item{\dots}{Other arguments to pass to \code{\link{gmmFit}} or
    \code{\link{gelFit}}.}
}

\section{Methods}{
\describe{

\item{\code{signature(object = "ANY")}}{
The method from the \pkg{stats} in used in that case.
}

\item{\code{signature(object = "gelfit")}}{
Method for any GEL fit class.
}

\item{\code{signature(object = "gmmfit")}}{
Method for any GMM fit class.
}

\item{\code{signature(object = "numeric")}}{
It computes the GEL confidence interval for the mean.
}

\item{\code{signature(object = "data.frame")}}{
  It computes the 2D GEL confidence region for the means of two
  variables.
}

\item{\code{signature(object = "matrix")}}{
  It converts the object into a data.frame and call its method.
}

}}
\keyword{methods}
\keyword{Confidence Intervals}
