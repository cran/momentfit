## ----echo=FALSE---------------------------------------------------------------
library(knitr)
opts_chunk$set(size='footnotesize')

## -----------------------------------------------------------------------------
library(momentfit)
data(simData)
lin <- momentModel(y~x1+x2, ~x2+z1+z2, data=simData, vcov="iid")  

## -----------------------------------------------------------------------------
theta0=c(mu=1,sig=1)
x <- simData$x1
dat <- data.frame(x=x, x2=x^2, x3=x^3, x4=x^4)
gform <- list(x~mu,
              x2~mu^2+sig, 
              x3~mu^3+3*mu*sig,
              x4~mu^4+6*mu^2*sig+3*sig^2)
form <- momentModel(gform, NULL, theta0=theta0, vcov="MDS", data=dat)
form

## -----------------------------------------------------------------------------
fct <- function(theta, x)
    cbind(x-theta[1], (x-theta[1])^2-theta[2],
          (x-theta[1])^3, (x-theta[1])^4-3*theta[2]^2)
dfct <- function(theta, x)
    {
        m1 <- mean(x-theta[1])
        m2 <- mean((x-theta[1])^2)
        m3 <- mean((x-theta[1])^3)
        matrix(c(-1, -2*m1, -3*m2, -4*m3, 
                 0, -1, 0, -6*theta[2]), 4, 2)
    }
theta0=c(mu=1,sig2=1)
func <- momentModel(fct, simData$x3, theta0=theta0, grad=dfct, vcov="iid")
func

## -----------------------------------------------------------------------------
theta0=c(b0=1, b1=1, b2=1)
gform <- y~exp(b0+b1*x1+b2*x2)
nlin <- momentModel(gform, ~z1+z2+z3+x2, theta0=theta0, vcov="MDS", data=simData)
nlin

## -----------------------------------------------------------------------------
rhoEL

## -----------------------------------------------------------------------------
args(getLambda)

## -----------------------------------------------------------------------------
X <- simData[c("x3","x4","z5")]
(res <- getLambda(X, gelType="EL"))$lambda
res$convergence$convergence

## -----------------------------------------------------------------------------
(res <- getLambda(X, gelType="ET",control=list(maxit=2000)))$lambda
res$convergence$convergence

## -----------------------------------------------------------------------------
(res <- getLambda(X, rhoFct=rhoEEL))$lambda
res$convergence$convergence

## -----------------------------------------------------------------------------
showMethods("solveGel")

## -----------------------------------------------------------------------------
mylSolve <- function(gmat, lambda0, gelType=NULL, rhoFct=NULL, k=1, ...) 
    {
        lambda <-  rep(0,ncol(gmat))
        obj <-  sum(colMeans(gmat)^2)
        list(lambda=lambda, convergence=0, obj=obj)
    }

## -----------------------------------------------------------------------------
solveGel(lin,theta0=c(0,0,0), lamSlv=mylSolve)

## -----------------------------------------------------------------------------
mylSolve <- function(gmat, lambda0=NULL, gelType=NULL, rhoFct=NULL, k=1, ...) 
    {
        gmat <- as.matrix(gmat)
        res  <-  getLambda(gmat, lambda0, gelType="ET", k=k)
        gml <- c(gmat%*%res$lambda)
        w <- exp(gml)
        w <- w/sum(w)
        n <- nrow(gmat)
        res$obj <- mean(-log(w*n))
        res
    }
etelFit <- solveGel(lin,theta0=c(1,1,0), lamSlv=mylSolve)
etelFit$theta

## -----------------------------------------------------------------------------
solveGel(update(lin, gelType="ETEL"), theta0=c(1,1,0))$theta

## ----warning=FALSE------------------------------------------------------------
solveGel(lin, theta0=c(1,1,0), lControl=list(algo="nlminb"))$theta

## -----------------------------------------------------------------------------
res <- solveGel(lin, theta0=c(1,1,0), lControl=list(restrictedLam=c(2L,3L)))
res$lambda

## -----------------------------------------------------------------------------
solveGel(lin, theta0=c(1,1,0),
         tControl=list(method="BFGS", control=list(maxit=2000, reltol=1e-9)))$theta

## ----warning=FALSE------------------------------------------------------------
fit <- gelFit(lin)

## -----------------------------------------------------------------------------
showClass("gelfit")

## -----------------------------------------------------------------------------
fit <- gelFit(lin)
print(fit, lambda=FALSE, model=FALSE)

## -----------------------------------------------------------------------------
pt <- getImpProb(fit)
pt$convMom

## -----------------------------------------------------------------------------
summary(fit)

## -----------------------------------------------------------------------------
specTest(fit)

## -----------------------------------------------------------------------------
confint(fit, parm=1:2, level=0.90)
confint(fit, level=0.90, lambda=TRUE)

## -----------------------------------------------------------------------------
confint(fit, 1:2, type="invLR", level=0.90)

## -----------------------------------------------------------------------------
muModel <- momentModel(x1~1, ~1, data=simData)

## -----------------------------------------------------------------------------
muFit <- gelFit(muModel, gelType="EL", tControl=list(method="Brent", lower=-10, upper=10))
muFit@theta

## -----------------------------------------------------------------------------
mean(simData$x1)

## -----------------------------------------------------------------------------
confint(muFit, type="invLR")

## -----------------------------------------------------------------------------
fit <- gelFit(lin)
fit@theta

## -----------------------------------------------------------------------------
fit2 <- update(fit, coefSlv="nlminb")
fit2@theta

## -----------------------------------------------------------------------------
fit3 <- update(fit2, gelType="ET", lControl=list())
fit3@theta

## -----------------------------------------------------------------------------
rlin <- restModel(fit3@model, "x1=1")
update(fit3, newModel=rlin)

## -----------------------------------------------------------------------------
mod1 <- momentModel(y~x1+x2, ~x2+z1+z2+z3, data=simData)
fit1 <- gelFit(mod1)

## -----------------------------------------------------------------------------
eta <- c(coef(fit1), fit1@lambda)
names(eta) <- NULL

## -----------------------------------------------------------------------------
mod2 <-  momentModel(momFct, fit1, theta0=eta, vcov="MDS")

## -----------------------------------------------------------------------------
fit2 <- evalGmm(mod2, eta)
v <- vcov(fit2)[1:3,1:3]
sqrt(diag(v))

## -----------------------------------------------------------------------------
sqrt(diag(vcov(fit1)$vcov_par))

## -----------------------------------------------------------------------------
sqrt(diag(vcov(fit1, robToMiss=TRUE)$vcov_par))

## -----------------------------------------------------------------------------
summary(fit1, robToMiss=TRUE)@coef

## -----------------------------------------------------------------------------
print(evalGel(lin, theta=c(1,2,3), lambda=rep(.1,4)), model=FALSE)

## -----------------------------------------------------------------------------
specTest(evalGel(muModel, theta=4), type="LR")

## -----------------------------------------------------------------------------
rmuModel <- restModel(muModel, R=matrix(1,1,1), rhs=4)
specTest(gelFit(rmuModel))

## -----------------------------------------------------------------------------
fit <- gel4(y~x1+x2, ~x2+z1+z2, data=simData, vcov="iid", gelType="EL",
             lControl=list(algo="Wu"))  
print(fit, model=FALSE)

## -----------------------------------------------------------------------------
fit <- gel4(y~x1+x2, ~x2+z1+z2, data=simData, vcov="iid", gelType="EL",
             lControl=list(algo="Wu"), theta0=c(1,1,0), initTheta="theta0")  
coef(fit)

## -----------------------------------------------------------------------------
fit <- gel4(y~x1+x2, ~x2+z1+z2, data=simData, vcov="iid", gelType="EL",
             lControl=list(algo="Wu"), cstLHS="x2=x1")
print(fit, model=FALSE)

## -----------------------------------------------------------------------------
theta0=c(mu=1,sig=1)
x <- simData$x1
dat <- data.frame(x=x, x2=x^2, x3=x^3, x4=x^4)
gform <- list(x~mu,
              x2~mu^2+sig, 
              x3~mu^3+3*mu*sig,
              x4~mu^4+6*mu^2*sig+3*sig^2)
fit <- gel4(g=gform, gelType="EEL", theta0=theta0, vcov="MDS", data=dat)
print(fit, model=FALSE)

## -----------------------------------------------------------------------------
fct <- function(theta, x) 
   cbind(x-theta[1], (x-theta[1])^2-theta[2],
          (x-theta[1])^3, (x-theta[1])^4-3*theta[2]^2)
dfct <- function(theta, x)
    {
        m1 <- mean(x-theta[1])
        m2 <- mean((x-theta[1])^2)
        m3 <- mean((x-theta[1])^3)
        matrix(c(-1, -2*m1, -3*m2, -4*m3, 
                 0, -1, 0, -6*theta[2]), 4, 2)
    }
fit <- gel4(g=fct, x=simData$x3, theta0=c(1,1), grad=dfct, vcov="iid", gelType="ET")
print(fit, model=FALSE)

## -----------------------------------------------------------------------------
theta0=c(b0=1, b1=0, b2=0)
gform <- y~exp(b0+b1*x1+b2*x2)
fit <- gel4(gform, ~z1+z2+z3+x2, theta0=theta0, vcov="MDS", data=simData,
                 gelType="HD")
print(fit, model=FALSE)

## -----------------------------------------------------------------------------
x2 <- simData$x2
confint(x2, gelType="EL")
print(confint(x2, gelType="EL", type="invLR"), digits=5)

## -----------------------------------------------------------------------------
confint(simData, "x2", gelType="EL")

## -----------------------------------------------------------------------------
res <- confint(simData, c("x2","y"), npoints=20)
res

## ----fig.height=5, eval=FALSE-------------------------------------------------
#  res2 <- confint(simData, c("x2","y"), type="invLR", npoints=20)
#  c1 <- col2rgb("darkorange2")/255
#  c1 <- rgb(c1[1],c1[2],c1[3],.5)
#  c2 <- col2rgb("lightblue2")/255
#  c2 <- rgb(c2[1],c2[2],c2[3],.5)
#  plot(res, pch=20, bg=1, Pcol=c1, col=c1, density=10, ylim=c(3.5,6.5),
#       xlim=c(4.8,7.5))
#  plot(res2, pch=20, bg=2, Pcol=c2, col=c2, density=10, add=TRUE)
#  legend("topright", c("Wald","LR"), col=c(c1,c2), pch=20 , bty="n")

## ----extract, message=FALSE, warning=FALSE------------------------------------
library(texreg)
setMethod("extract", "gelfit", 
          function(model, includeSpecTest=TRUE, 
                   specTest=c("LR","LM","J"), include.nobs=TRUE, 
                   include.obj.fcn=TRUE, ...)
              {
                  specTest <- match.arg(specTest)
                  s <- summary(model, ...)
                  wspecTest <- grep(specTest, rownames(s@specTest@test))
                  spec <- modelDims(model@model)
                  coefs <- s@coef
                  names <- rownames(coefs)
                  coef <- coefs[, 1]
                  se <- coefs[, 2]
                  pval <- coefs[, 4]
                  n <- model@model@n
                  gof <- numeric()
                  gof.names <- character()
                  gof.decimal <- logical()
                  if (includeSpecTest) {
                      if (spec$k == spec$q)
                          {
                              obj.fcn <- NA
                              obj.pv <- NA
                          } else {
                              obj.fcn <- s@specTest@test[wspecTest,1]
                              obj.pv <- s@specTest@test[wspecTest,3]
                          }
                      gof <- c(gof, obj.fcn, obj.pv)                      
                      gof.names <- c(gof.names, 
                                     paste(specTest,"-test Statistics", sep=""),
                                     paste(specTest,"-test p-value", sep=""))
                      gof.decimal <- c(gof.decimal, TRUE, TRUE)
                  }
                  if (include.nobs == TRUE) {
                      gof <- c(gof, n)
                      gof.names <- c(gof.names, "Num.\\ obs.")
                      gof.decimal <- c(gof.decimal, FALSE)
                  }
                  tr <- createTexreg(coef.names = names, coef = coef, se = se, 
                                     pvalues = pval, gof.names = gof.names, gof = gof, 
                                     gof.decimal = gof.decimal)
                  return(tr)
              })

## ----results='asis'-----------------------------------------------------------
fit1 <- gel4(y~x1, ~x2+x3+x4, data=simData)
fit2 <- gel4(y~x1+x2, ~x2+x3+x4, data=simData)
texreg(list(fit1,fit2), digits=4)

