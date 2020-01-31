\name{systemGmm-doc}
\alias{systemGmm}

\title{A guide to estimating systems of equations}
\description{
This document is meant to describe how to create system of equations
objects, estimating them and peforming hypothesis tests.
}
\details{
Instread of repeating the same example for each method, we are going
through all methods and classes for systems of equations.
}
  
\examples{
data(simData)

## first, we create an sysGmm object
g1 <- y1~x1+x4; h1 <- ~x4+z1+z2+z3+z4
g2 <- y2~x1+x2+x3; h2 <- ~x3+z1+z2+z3+z4
g3 <- y3~x2+x3+x4; h3 <- ~x3+x4+z1+z2+z3+z4
g <- list(g1,g2,g3)
h <- list(h1,h2,h3)
smodel <- sysMomentModel(g, h, data=simData, vcov="MDS")

## The show or print method
smodel

## The ']' method
smodel[1:2]
smodel[1] ## becomes a one equation model

## equation by equation 2SLS
tsls(smodel)

## or manually
lapply(1:3, function(i) coef(tsls(smodel[i])))

## Fitting the model by two-step GMM
res <- gmmFit(smodel)

## testing Overidentifying restrictions
specTest(res)

## All info using the summary method
## which includes equation by equation measures of
## the instrument stengths
summary(res)

### When the error id iid (homoscedastic), we have a
### FIVE estimator with 2SLS  as the first step
smodel <- sysMomentModel(g, h, data=simData, vcov="iid")
gmmFit(smodel)

### When the error is iid (homoscedastic), 
### all instruments are the same, and the first step is 2SLS,
### we have 3SLS
smodel <- sysMomentModel(g, ~x4+z1+z2+z3+z4, data=simData, vcov="iid")
gmmFit(smodel, initW='tsls')

### When the error is iid (homoscedastic), 
### the instruments are the same and are the union of all regressors,
### we have SUR
smodel <- sysMomentModel(g, NULL, data=simData, vcov="iid")
gmmFit(smodel, initW='tsls')

############ Restricted models ##################

## unrestricted
smodel <- sysMomentModel(g, h, data=simData, vcov="MDS")
res <- gmmFit(smodel)

## no cross-equation restrictions
R1 <- list(c("x1=-12*x4"), character(), c("x2=0.8", "x4=0.3"))
rm1 <- restModel(smodel, R1)
(res1 <- gmmFit(rm1))

## Cross equation restrictions
R2<- c("Eqn1.x1=1", "Eqn2.x1=Eqn3.x2")
rm2 <- restModel(smodel, R2)
(es2 <- gmmFit(rm2))## no longer expressed as a system

## testing the restriction

hypothesisTest(res, res1, type="LR")
hypothesisTest(res, res1, type="LM")
hypothesisTest(res, res1, type="Wald")

}
\keyword{FIVE}
\keyword{3SLS}
\keyword{SUR}
\keyword{GMM}
\keyword{System of Equations}
