\name{Mroz}
\alias{Mroz}
\docType{data}
\title{
Labour data on married women
}
\description{
The dataset was used by Mroz (1987) and in examples in Wooldridge (2016)
}
\usage{data("Mroz")}
\format{
  A data frame with 753 observations on the following 22 variables.
  \describe{
    \item{\code{inlf}}{=1 if in lab frce, 1975}
    \item{\code{hours}}{hours worked, 1975}
    \item{\code{kidslt6}}{number of kids < 6 years}
    \item{\code{kidsge6}}{number of kids 6-18}
    \item{\code{age}}{woman's age in years}
    \item{\code{educ}}{years of schooling}
    \item{\code{wage}}{Estimated wage from earnings and hours}
    \item{\code{repwage}}{reported wage at interview in 1976}
    \item{\code{hushrs}}{hours worked by husband, 1975}
    \item{\code{husage}}{husband's age}
    \item{\code{huseduc}}{husband's years of schooling}
    \item{\code{huswage}}{husband's hourly wage, 1975}
    \item{\code{faminc}}{family income, 1975}
    \item{\code{mtr}}{federal marginal tax rate facing woman}
    \item{\code{motheduc}}{mother's years of schooling}
    \item{\code{fatheduc}}{father's years of schooling}
    \item{\code{unem}}{unemployment rate in county of residence}
    \item{\code{city}}{=1 if live in SMSA}
    \item{\code{exper}}{actual labor market experience}
    \item{\code{nwifeinc}}{(\eqn{faminc - wage*hours)/1000}}
  }
}

\source{
From Wooldridge (2016) online resources.
}
\references{
Mroz, T.A. (1987), \emph{The Sensitivity of an Empirical Model of
  Married Women's Hours of Work to Economic and Statistical Assumptions},
 Econometrica, \bold{55}, 657--678.
387--405.

Wooldridge, J.M. (2016). \emph{Introductory Econometrics, A Modern
  Approach, 6th edition}, Cengage Learning.  
}
\examples{

## Example 15.1 of Wooldridge (2016)

data(Mroz)
Mroz <- subset(Mroz, hours>0)
## I guess IID is assumed (That's how we get the same s.e.)
## By default a sandwich vcov is computed because it is 
## a just-identified model.
res4 <- gmm4(log(wage)~educ, ~fatheduc, vcov="iid", data=Mroz)
summary(res4)

## If we adjust the variance of the residuals, however,
## we are a little off (very little)

summary(res4, df.adj=TRUE)


## Example 15.5 of Wooldridge (2016)
## Need to adjust for degrees of freedom in order
## to get the same s.e.
## The first stage F-test is very different though
## Cannot get the same even if do it manually
## with the linearHypothesis from the car package
model <- momentModel(log(wage)~educ+exper+I(exper^2),
~exper+I(exper^2)+fatheduc+motheduc, vcov="iid", data=Mroz)
res <- tsls(model)
summary(res, df.adj=TRUE)

}
\keyword{datasets}
