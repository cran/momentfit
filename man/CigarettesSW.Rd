\name{CigarettesSW}
\alias{CigarettesSW}
\title{Cigarette Consumption Panel Data}

\description{
Panel data on cigarette consumption for the 48 continental US States from 1985--1995.
}

\usage{data("CigarettesSW")}

\format{
A data frame containing 48 observations on 7 variables for 2 periods. 
\describe{
\item{state}{Factor indicating state.}
\item{year}{Factor indicating year.}
\item{cpi}{Consumer price index.}
\item{population}{State population.}
\item{packs}{Number of packs per capita.}
\item{income}{State personal income (total, nominal).}
\item{tax}{Average state, federal and average local excise taxes for fiscal year.}
\item{price}{Average price during fiscal year, including sales tax.}
\item{taxs}{Average excise taxes for fiscal year, including sales tax.}
}
}

\source{
Online complements to Stock and Watson (2007). The dataset and this help file comes from
the AER package.
}

\references{
Stock, J.H. and Watson, M.W. (2007). \emph{Introduction to Econometrics},
2nd ed. Boston: Addison Wesley.

Christian Kleiber and Achim Zeileis (2008). Applied Econometrics with
R. New York: Springer-Verlag. ISBN 978-0-387-77316-2. URL
https://CRAN.R-project.org/package=AER
}

\examples{
## Stock and Watson (2007)
## data and transformations
data(CigarettesSW)
CigarettesSW$rprice <- with(CigarettesSW, price/cpi)
CigarettesSW$rincome <- with(CigarettesSW, income/population/cpi)
CigarettesSW$tdiff <- with(CigarettesSW, (taxs - tax)/cpi)
c1985 <- subset(CigarettesSW, year == "1985")
c1995 <- subset(CigarettesSW, year == "1995")

## Equation 12.15
model1 <- momentModel(log(packs)~log(rprice)+log(rincome),
                   ~log(rincome)+tdiff, data = c1995, vcov="MDS")
res1 <- gmmFit(model1)

## HC0 robust se (different from the textbook)
summary(res1, sandwich=TRUE)

## HC1 robust se (like in the textbook)
## A little harder to get, but is it really worth it
## in the case of GMM?

summary(res1, sandwich=TRUE, df.adj=TRUE)@coef

## Equation 12.16
model2<- momentModel(log(packs)~log(rprice)+log(rincome),
                  ~log(rincome)+tdiff+I(tax/cpi), data = c1995,
                  centeredVcov=FALSE, vcov="MDS")
res2<- tsls(model2)
summary(res2, sandwich=TRUE, df.adj=TRUE)

## Table 12.1
data <- data.frame(dQ=log(c1995$pack/c1985$pack),
                   dP=log(c1995$rprice/c1985$rprice),
                   dTs=c1995$tdiff-c1985$tdiff,
                   dT=c1995$tax/c1995$cpi-c1985$tax/c1985$cpi,
                   dInc=log(c1995$rincome/c1985$rincome))
model1 <- momentModel(dQ~dP+dInc, ~dInc+dTs, vcov="MDS", data=data)
model2 <- momentModel(dQ~dP+dInc, ~dInc+dT, vcov="MDS", data=data)
model3 <- momentModel(dQ~dP+dInc, ~dInc+dTs+dT, vcov="MDS", data=data)

res1 <- tsls(model1)
summary(res1, TRUE, TRUE)
res2 <- tsls(model2)
summary(res2, TRUE, TRUE)
res3 <- tsls(model3)
summary(res3, TRUE, TRUE)
}

\keyword{datasets}
