\name{LabourCR}
\alias{LabourCR}
\docType{data}
\title{
Labour data from Greene (2012) applications,
}
\description{
A panel data set of 565 individuals from 1976 to 1982 used by Cornwell
and Rupert (1988)
}
\usage{data("LabourCR")}
\format{
  A data frame with 4165 observations on the following 12 variables.
  \describe{
    \item{\code{EXP}}{Year of full time experience.}
    \item{\code{WKS}}{Weeks worked.}
    \item{\code{OCC}}{1 if blue-collar occupation, 0 otherwise.}
    \item{\code{IND}}{1 if works in a manufacture industry, 0 otherwise.}
    \item{\code{SOUTH}}{1 if resides in the south, 0 otherwise.}
    \item{\code{SMSA}}{1 if resides in an SMSA, 0 otherwise.}
    \item{\code{MS}}{1 if married, 0 otherwise.}
    \item{\code{FEM}}{1 if the individual is a female and 0 otherwise.}
    \item{\code{UNION}}{1 if wage is set by a union contract and 0 otherwise.}
    \item{\code{ED}}{Years of education.}
    \item{\code{BLK}}{1 if the individual is black and 0 otherwise.}
    \item{\code{LWAGE}}{Log wage.}
  }
}

\source{
  Greene (2012) online resources:
  (http://pages.stern.nyu.edu/~wgreene/Text/Edition7/tablelist8new.htm)
}

\references{
Green, W.H.. (2012). \emph{Econometric Analysis, 7th edition}, Prentice Hall.  

Cornwell, C. and Rupert, P. (1988), \emph{Efficient Estimation with
  Panel Data: An Empirical Comparision of Instrumental Variable
  Estimators}, Journal of Applied Econometrics, No.3, 149--155.
}
\examples{
data(LabourCR)
## Table 8.1 of Greene (2012)
## Model with Z2 (iid is assumed in Table 8.1 given the s.e.)
model2 <- momentModel(WKS~LWAGE+ED+UNION+FEM, ~IND+ED+UNION+FEM+SMSA, vcov="iid",
                   data=LabourCR)
## Model with Z1 using the subsetting method '['
model1 <- model2[-6L]

# Second column
res1 <- tsls(model1)
summary(res1)@coef

# Third column
res2 <- tsls(model2)
summary(res2)@coef

}
\keyword{datasets}
