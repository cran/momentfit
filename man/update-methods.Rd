\name{update-methods}
\docType{methods}
\alias{update-methods}
\alias{update,ANY-method}
\alias{update,list-method}
\alias{update,gmmfit-method}
\alias{update,momentModel-method}
\alias{update,gelfit-method}
\title{ ~~ Methods for Function \code{update} in Package \pkg{stats} ~~}
\description{
The method is used to refit a model with either a different
method or with modifications to the \code{momentModel}.
}
\usage{
\S4method{update}{gmmfit}(object, \dots, evaluate=TRUE)

\S4method{update}{momentModel}(object, \dots)

\S4method{update}{gelfit}(object, newModel=NULL, \dots,
evaluate=TRUE)

\S4method{update}{list}(object, \dots)

}
\arguments{
  \item{object}{An object produced by \code{"gelFit"}, \code{"gmmFit"}
    or a model. It can also be a list, in which case, it is used to
    change elements of a list. }  
  \item{\dots}{Arguments to modify the model or the GMM method}
  \item{newModel}{When provided, the new model is estimated using the
    same specification. For example, it is particularly useful to
    estimate the restricted model using the same optim specification as
    the unrestricted model.}
  \item{evaluate}{The modified \code{call} argument is only evaluated
    when evaluate is \code{TRUE}}
  }
\section{Methods}{
\describe{

  \item{\code{signature(object = "ANY")}}{
    That just calls \code{"update"} from the \code{"stats"} package.
}

\item{\code{signature(object = "gmmfit")}}{
}

\item{\code{signature(object = "momentModel")}}{
}

\item{\code{signature(object = "list")}}{
}
}}
\examples{
x <- rchisq(200,5)
z1 <- rnorm(200)
z2 <- .2*x+rnorm(200)
y <- x+rnorm(200)
dat <- data.frame(y=y,z1=z1,x=x,z2=z2)
theta <- c(beta0=1,beta1=2)
model1 <- momentModel(y~x, ~z1+z2, data=dat)

(res <- gmmFit(model1))

## lets change to iterative
update(res, type="iter")

## Let change the HAC specification in the model1 object
## to MDS
update(res, vcov="MDS")


}
\keyword{methods}
