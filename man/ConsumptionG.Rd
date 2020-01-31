\name{ConsumptionG}
\alias{ConsumptionG}
\docType{data}
\title{
Consumption data from Greene (2012) applications.
}
\description{
Quarterly macroeconomic US data from 1950 to 2000.
}
\usage{data("ConsumptionG")}
\format{
  A data frame with 204 observations on the following 14 variables.
  \describe{
    \item{\code{YEAR}}{Year}
    \item{\code{QTR}}{Quarter}
    \item{\code{REALGDP}}{Read GDP}
    \item{\code{REALCONS}}{Real Consumption}
    \item{\code{REALINVS}}{Real Investment}
    \item{\code{REALGOVT}}{Real public expenditure}
    \item{\code{REALDPI}}{ector}
    \item{\code{CPI_U}}{CPI}
    \item{\code{M1}}{Money stock}
    \item{\code{TBILRATE}}{Interest rate}
    \item{\code{UNEMP}}{Unemployment rate}
    \item{\code{POP}}{Population}
    \item{\code{INFL}}{Inflation}
    \item{\code{REALINT}}{Real interest rate.}
  }
}

\source{
  Greene (2012) online resources:
  (http://pages.stern.nyu.edu/~wgreene/Text/Edition7/tablelist8new.htm)  
}
\references{
Green, W.H.. (2012). \emph{Econometric Analysis, 7th edition}, Prentice Hall.    
}
\examples{
data(ConsumptionG)
## Get the data ready for Table 8.2 of Greene (2012)
Y <- ConsumptionG$REALDPI
C <- ConsumptionG$REALCONS
n <- nrow(ConsumptionG)
Y1 <- Y[-c(1,n)]; Y2 <- Y[-c(n-1,n)]; Y <- Y[-c(1:2)]
C1 <- C[-c(1,n)]; C <- C[-(1:2)]
dat <- data.frame(Y=Y,Y1=Y1,Y2=Y2,C=C,C1=C1)

## Starting at the NLS estimates (from the table)
theta0=c(alpha=468, beta=0.0971, gamma=1.24)

## Greene (2012) seems to assume iid errors (probably wrong assumption here)
model <- momentModel(C~alpha+beta*Y^gamma, ~C1+Y1+Y2, data=dat, theta0=theta0, vcov="iid")

### Scaling the parameters increase the speed of convergence
res <- gmmFit(model, control=list(parscale=c(1000,.1,1)))

### It also seems that there is a degree of freedom adjustment for the
### estimate of the variance of the error term.
summary(res, df.adj=TRUE)@coef



}
\keyword{datasets}
