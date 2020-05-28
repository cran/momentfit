## ----echo=FALSE---------------------------------------------------------------
library(knitr)
opts_chunk$set(size='footnotesize')

## ----warning=FALSE, message=FALSE---------------------------------------------
library(momentfit)
data(simData)
modelF <- model.frame(y~x1+x2, simData)
instF <- model.frame(~x2+z1+z2, simData)
mod1 <- new("linearModel", modelF=modelF, instF=instF, k=3L, q=4L, vcov="iid",
             parNames=c("(Intercept)", "x1","x2"), n=50L, 
             momNames=c("(Intercept)", "x2", "z1", "z2"), 
             isEndo=c(FALSE, TRUE, FALSE, FALSE), smooth=FALSE)

## -----------------------------------------------------------------------------
mod1

## -----------------------------------------------------------------------------
mod1 <- momentModel(y~x1+x2, ~x2+z1+z2, data=simData, vcov="iid")
mod1

## -----------------------------------------------------------------------------
theta0 <- c(theta0=1, theta1=1, theta2=2)
mod2 <- momentModel(y~exp(theta0+theta1*x1+theta2*x2), ~x2+z1+z2, theta0, 
                 data=simData, vcov="iid")
mod2

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

## -----------------------------------------------------------------------------
theta0=c(mu=1,sig2=1)
x <- simData$x3
mod3 <- momentModel(fct, x, theta0, grad=dfct, vcov="iid")
mod3

## -----------------------------------------------------------------------------
theta0=c(mu=1,sig=1)
dat <- data.frame(x=x, x2=x^2, x3=x^3, x4=x^4)
gform <- list(x~mu,
              x2~mu^2+sig, 
              x3~mu^3+3*mu*sig,
              x4~mu^4+6*mu^2*sig+3*sig^2)
mod4 <- momentModel(gform, NULL, theta0, vcov="MDS", data=dat)
mod4

## ----eval=FALSE---------------------------------------------------------------
#  dat <- data.frame(x=x)
#  gform <- list(x~mu,
#                x^2~mu^2+sig,
#                x^3~mu^3+3*mu*sig,
#                x^4~mu^4+6*mu^2*sig+3*sig^2)
#  mod4 <- momentModel(gform, NULL, theta0, vcov="MDS", data=dat)

## -----------------------------------------------------------------------------
mod.hac <- momentModel(y~x1+x2, ~x1+z2+z3, vcov="HAC", 
                    vcovOptions=list(kernel="Bartlett", bw="NeweyWest"),
                    data=simData)
mod.hac

## -----------------------------------------------------------------------------
data("InstInnovation", package = "sandwich")

## -----------------------------------------------------------------------------
mod.cl1 <- momentModel(sales~value, ~value, vcov="CL", 
                    vcovOptions=list(cluster=~company),
                    data=InstInnovation)
mod.cl1

## -----------------------------------------------------------------------------
mod.cl2 <- momentModel(sales~value, ~value, vcov="CL", 
                    vcovOptions=list(cluster=~company+year, multi0=TRUE),
                    data=InstInnovation)
mod.cl2

## -----------------------------------------------------------------------------
smod1 <- momentModel(y~x1+x2, ~x2+z1+z2, data=simData, smooth=TRUE)  
smod1

## -----------------------------------------------------------------------------
smod2 <- momentModel(y~x1+x2, ~x2+z1+z2, data=simData, smooth=TRUE,
                     vcovOptions=list(kernel="Parzen", bw="NeweyWest", prewhite=1))  
smod2

## -----------------------------------------------------------------------------
smod2@sSpec

## -----------------------------------------------------------------------------
smod2@sSpec@w

## -----------------------------------------------------------------------------
setCoef(mod2, 1:3)

## -----------------------------------------------------------------------------
setCoef(mod2, c(theta1=1, theta2=1, theta0=2))

## -----------------------------------------------------------------------------
theta0 <- c(theta0=1, theta1=1, theta2=2)
e1 <- residuals(mod1, c(1,2,3))
e2 <- residuals(mod2, theta0)

## -----------------------------------------------------------------------------
e1 <- Dresiduals(mod1)
theta0 <- setCoef(mod2, c(1,1,2))
e2 <- Dresiduals(mod2, theta0)

## -----------------------------------------------------------------------------
Z <- model.matrix(mod1, type="instruments")

## -----------------------------------------------------------------------------
X <- model.matrix(mod1)

## -----------------------------------------------------------------------------
Y <- modelResponse(mod1)

## -----------------------------------------------------------------------------
mod1[1:3]
mod2[c(1,2,4)]
mod3[-1]

## -----------------------------------------------------------------------------
mod4 <- as(mod1, "nonlinearModel")

## -----------------------------------------------------------------------------
mod4@parNames
mod4@fLHS
mod4@fRHS

## -----------------------------------------------------------------------------
subset(mod1, simData$x1>4)

## -----------------------------------------------------------------------------
gt <- evalMoment(mod1, 1:3)

## -----------------------------------------------------------------------------
theta0 <- c(theta0=.1, theta1=1, theta2=-2)
## or ##
theta0 <- setCoef(mod2, c(.1,1,-2))
evalDMoment(mod2, theta0)

## -----------------------------------------------------------------------------
vcov(mod1, theta=c(1,1,1))

## -----------------------------------------------------------------------------
momentStrength(mod1)

## -----------------------------------------------------------------------------
update(mod1, vcov="MDS")

## -----------------------------------------------------------------------------
update(mod1, vcov="HAC")

## -----------------------------------------------------------------------------
update(mod1, smooth=TRUE)

## -----------------------------------------------------------------------------
UR.mod1 <- momentModel(y~x1+x2+x3+z1, ~x1+x2+z1+z2+z3+z4, data=simData)

## -----------------------------------------------------------------------------
R1 <- matrix(c(1,1,0,0,0,0,0,2,0,0,0,0,0,1,-1),3,5, byrow=TRUE)
q1 <- c(0,1,3)
R1.mod1 <- restModel(UR.mod1, R1, q1)
R1.mod1

## -----------------------------------------------------------------------------
R2 <- c("x1","2*x2+z1=2", "4+x3*5=3")
R2.mod1 <- restModel(UR.mod1, R2)
printRestrict(R2.mod1)

## -----------------------------------------------------------------------------
UR.mod2 <- momentModel(y~x1*x2+exp(x3)+I(z1^2), ~x1+x2+z1+z2+z3+z4, data=simData)
R3 <- c("x1","exp(x3)+2*x1:x2", "I(z1^2)=3")
R3.mod2 <- restModel(UR.mod2, R3)
printRestrict(R3.mod2)

## -----------------------------------------------------------------------------
R1 <- c("theta1=theta2^2")
restModel(mod2, R1)
printRestrict(restModel(mod2, theta1~theta2))

## -----------------------------------------------------------------------------
restModel(mod3, "mu=0.5")

## -----------------------------------------------------------------------------
printRestrict(R2.mod1)

## -----------------------------------------------------------------------------
coef(R2.mod1, c(1.5,.5))

## -----------------------------------------------------------------------------
setCoef(R2.mod1, c(1.5,.5))

## -----------------------------------------------------------------------------
e1 <- residuals(as(R2.mod1, "linearModel"), 
               coef(R2.mod1, c(1.5,.5)))

## -----------------------------------------------------------------------------
e2 <- residuals(R2.mod1, c(1.5,.5))
all.equal(e1,e2)

## -----------------------------------------------------------------------------
R1 <- c("theta1=theta2^2")
R1.mod2 <- restModel(mod2, R1)
evalDMoment(mod2, c(theta0=1, theta1=1, theta2=1))
## with setCoef:
evalDMoment(R1.mod2, setCoef(R1.mod2, c(1,1)))

## -----------------------------------------------------------------------------
mod2@parNames
R1.mod2@parNames

## -----------------------------------------------------------------------------
modelDims(mod2)$parNames
modelDims(mod2)$k
modelDims(R1.mod2)$parNames
modelDims(R1.mod2)$k

## -----------------------------------------------------------------------------
model <- momentModel(y~x1, ~z1+z2, data=simData, vcov="iid") ## lets create a simple model
wObj <- evalWeights(model, w="ident")

## -----------------------------------------------------------------------------
wObj

## -----------------------------------------------------------------------------
evalWeights(model, w=diag(3))

## -----------------------------------------------------------------------------
wObj <- evalWeights(model, theta=c(1,2))

## -----------------------------------------------------------------------------
wObj@type

## -----------------------------------------------------------------------------
model2 <- momentModel(y~x1, ~z1+z2, data=simData, vcov="HAC")
evalWeights(model2, c(1,2))@type

## -----------------------------------------------------------------------------
evalWeights(model, w=diag(3))@type

## -----------------------------------------------------------------------------
wObj <- evalWeights(model, theta=1:2)

## -----------------------------------------------------------------------------
G <- evalDMoment(model, theta=1:2)
gbar <- colMeans(evalMoment(model, theta=1:2))

## -----------------------------------------------------------------------------
quadra(wObj, gbar)

## -----------------------------------------------------------------------------
quadra(wObj, G, gbar)

## -----------------------------------------------------------------------------
quadra(wObj)

## -----------------------------------------------------------------------------
wObj[1:2]

## -----------------------------------------------------------------------------
theta0 <- 1:2
wObj <- evalWeights(model, theta0)
theta1 <- 3:4
evalGmmObj(model, theta1, wObj)

## -----------------------------------------------------------------------------
mod <- momentModel(y~x1, ~z1+z2, data=simData, vcov="MDS")

## -----------------------------------------------------------------------------
wObj0 <- evalWeights(mod, w="ident")
res0 <- solveGmm(mod, wObj0)
res0$theta

## -----------------------------------------------------------------------------
wObj1 <- evalWeights(mod, res0$theta)
res1 <- solveGmm(mod, wObj1)
res1$theta

## -----------------------------------------------------------------------------
solveGmm(as(mod, "nonlinearModel"), wObj1)$theta
solveGmm(as(mod, "functionModel"), wObj1)$theta

## -----------------------------------------------------------------------------
theta0 <- c(theta0=0, theta1=0, theta2=0)
mod2 <- momentModel(y~exp(theta0+theta1*x1+theta2*x2), ~x2+z1+z2, theta0, 
                 data=simData, vcov="MDS")
wObj0 <- evalWeights(mod2, w="ident")
res1 <- solveGmm(mod2, wObj0, control=list(maxit=2000))
res1
solveGmm(mod2, wObj0, method="Nelder", control=list(maxit=2000))
solveGmm(mod2, wObj0, algo="nlminb", control=list(iter.max=2000))

## -----------------------------------------------------------------------------
R1 <- c("theta1=theta2^2")
rmod2 <- restModel(mod2, R1)
res2 <- solveGmm(rmod2, wObj0, control=list(maxit=2000))
res2

## -----------------------------------------------------------------------------
coef(rmod2, res2$theta)

## -----------------------------------------------------------------------------
mod <- momentModel(y~x1, ~z1+z2, data=simData, vcov="MDS")
gmmFit(mod, type="onestep")
print(gmmFit(mod, type="twostep"), model=FALSE)
print(gmmFit(mod, type="iter"), model=FALSE)

## -----------------------------------------------------------------------------
theta0 <- c(theta0=0, theta1=0, theta2=0)
mod2 <- momentModel(y~exp(theta0+theta1*x1+theta2*x2), ~x2+z1+z2, theta0, 
                 data=simData, vcov="MDS")
res1 <- gmmFit(mod2)
print(res1, model=FALSE)
theta0 <- c(theta0=0.5, theta1=0.5, theta2=-0.5)
res2 <- gmmFit(mod2, theta0=theta0, control=list(reltol=1e-8))
print(res2, model=FALSE)

## -----------------------------------------------------------------------------
theta0=c(mu=1,sig=1)
x <- rnorm(2000, 4, 5)
dat <- data.frame(x=x, x2=x^2, x3=x^3, x4=x^4)
gform <- list(x~mu,
              x2~mu^2+sig, 
              x3~mu^3+3*mu*sig,
              x4~mu^4+6*mu^2*sig+3*sig^2)
mod4 <- momentModel(gform, NULL, theta0, vcov="MDS", data=dat)
mod4@isMDE
print(gmmFit(mod4), model=FALSE)

## -----------------------------------------------------------------------------
print(gmmFit(mod4[1:2]), model=FALSE)

## -----------------------------------------------------------------------------
mod <- momentModel(y~x1, ~z1+z2+z3, data=simData, vcov="MDS")
res <- gmmFit(mod)
specTest(res)

## -----------------------------------------------------------------------------
specTest(res, 3)

## -----------------------------------------------------------------------------
summary(res)

## -----------------------------------------------------------------------------
summary(res, breadOnly=TRUE)@coef

## -----------------------------------------------------------------------------
mod <- momentModel(y~x1+x2+x3+z1, ~x1+x2+z1+z2+z3+z4, data=simData, vcov="iid")
res <- gmmFit(mod)

## -----------------------------------------------------------------------------
R <- c("x1=1", "x2=x3", "z1=-0.7")
rmod <- restModel(mod, R)
printRestrict(rmod)

## -----------------------------------------------------------------------------
hypothesisTest(object.u=res, R=R)

## -----------------------------------------------------------------------------
rmod@cstLHS
rmod@cstRHS

## -----------------------------------------------------------------------------
res.r <- gmmFit(rmod)

## -----------------------------------------------------------------------------
hypothesisTest(object.r=res.r)

## -----------------------------------------------------------------------------
hypothesisTest(object.r=res.r, object.u=res)

## ----eval=FALSE---------------------------------------------------------------
#  hypothesisTest(object.r=res.r, object.u=res, type="LM")
#  hypothesisTest(object.r=res.r, object.u=res, type="Wald")
#  hypothesisTest(object.r=res.r, object.u=res, type="LR")

## -----------------------------------------------------------------------------
coef(res.r)

## -----------------------------------------------------------------------------
e <- residuals(res)
e.r <- residuals(res.r)

## -----------------------------------------------------------------------------
mod <- momentModel(y~x1, ~z1+z2, data=simData, vcov="iid")
res1 <- gmmFit(mod)
res2 <- lm(y~x1, simData)
DWH(res1,res2)

## -----------------------------------------------------------------------------
DWH(res1)

## -----------------------------------------------------------------------------
confint(res1, level=0.99)

## -----------------------------------------------------------------------------
mod <- momentModel(y~x1+x2+z1, ~x1+z1+z2+z3, data=simData, vcov="iid")
res2 <- gmmFit(mod)
ci <- confint(res2, 2:3, area=TRUE)
ci

## ----fig.height=5-------------------------------------------------------------
plot(ci, col="lightblue", density=20, Pcol=2, bg=2)

## -----------------------------------------------------------------------------
res <- gmmFit(mod1)
res
update(res, vcov="MDS") ## changing only the model
update(res, vcov="MDS", type="iter")

## -----------------------------------------------------------------------------
mod <- momentModel(y~x1, ~z1+z2+z3, data=simData, vcov="MDS")
res <- tsls(mod)
summary(res)@coef

## -----------------------------------------------------------------------------
res1 <- gmm4(y~x1+x2, ~x2+z1+z2+z3, type="twostep", vcov="MDS", data=simData)
res1

## -----------------------------------------------------------------------------
res2 <- gmm4(y~x1+x2, ~x2+z1+z2+z3, type="iter", vcov="MDS", data=simData)

## -----------------------------------------------------------------------------
res1.r <- gmm4(y~x1+x2, ~x2+z1+z2+z3, type="twostep", vcov="MDS", 
               data=simData, cstLHS="x1=x2")
res1.r

## -----------------------------------------------------------------------------
hypothesisTest(res1, res1.r, type="LR")

## -----------------------------------------------------------------------------
res3 <- tsls(y~x1+x2, ~x2+z1+z2+z3, vcov="MDS", data=simData)
res3

## -----------------------------------------------------------------------------
res3 <- gmm4(y~theta0+exp(theta1*x1+theta2*x2), ~x2+z1+z2+z3+z4, vcov="iid",
             theta0=c(theta0=1, theta1=0, theta2=0), data=simData)
res3

## -----------------------------------------------------------------------------
update(res3, data=simData[1:45,])

## -----------------------------------------------------------------------------
update(res3, x = ~x2+z1+z2+z3, cstLHS="theta1=theta2")

## -----------------------------------------------------------------------------
data(CigarettesSW)
CigarettesSW$rprice <- with(CigarettesSW, price/cpi)
CigarettesSW$rincome <- with(CigarettesSW, income/population/cpi)
CigarettesSW$tdiff <- with(CigarettesSW, (taxs - tax)/cpi)
c1985 <- subset(CigarettesSW, year == "1985")
c1995 <- subset(CigarettesSW, year == "1995")

## -----------------------------------------------------------------------------
res1 <- gmm4(log(packs)~log(rprice)+log(rincome),
             ~log(rincome)+tdiff, data = c1995, vcov="MDS")
summary(res1, sandwich=TRUE, df.adj=TRUE)@coef

## -----------------------------------------------------------------------------
res2<- tsls(log(packs)~log(rprice)+log(rincome),
            ~log(rincome)+tdiff+I(tax/cpi), data = c1995,
            centeredVcov=FALSE, vcov="MDS")
summary(res2, sandwich=TRUE, df.adj=TRUE)@coef

## ----extract, echo=FALSE, message=FALSE, warning=FALSE------------------------
library(texreg)
setMethod("extract", "gmmfit", 
          function(model, includeJTest=TRUE, includeFTest=TRUE, ...)
              {
                  s <- summary(model, ...)
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
                  if (includeJTest) {
                      if (spec$k == spec$q)
                          {
                              obj.fcn <- NA
                              obj.pv <- NA
                          } else {
                              obj.fcn <- s@specTest@test[1]
                              obj.pv <- s@specTest@test[3]
                          }
                      gof <- c(gof, obj.fcn, obj.pv)
                      gof.names <- c(gof.names, "J-test Statistics", "J-test p-value")
                      gof.decimal <- c(gof.decimal, TRUE, TRUE)
                  }
                  if (includeFTest) {
                      str <- s@strength$strength
                      if (is.null(str))
                          {
                              gof <- c(gof, NA)
                              gof.names <- c(gof.names, "First Stage F-stats")
                              gof.decimal <- c(gof.decimal, TRUE)
                          } else {
                              for (i in 1:nrow(str))
                                  {
                                      gof <- c(gof, str[i,1])
                                      gofn <- paste("First Stage F-stats(",
                                                    rownames(str)[i], ")", sep="")
                                      gof.names <- c(gof.names, gofn)
                                      gof.decimal <- c(gof.decimal, TRUE)
                                  }
                          }
                  }
                  tr <- createTexreg(coef.names = names, coef = coef, se = se, 
                                     pvalues = pval, gof.names = gof.names, gof = gof, 
                                     gof.decimal = gof.decimal)
                  return(tr)
              })

## -----------------------------------------------------------------------------
data <- data.frame(dQ=log(c1995$pack/c1985$pack),
                   dP=log(c1995$rprice/c1985$rprice),
                   dTs=c1995$tdiff-c1985$tdiff,
                   dT=c1995$tax/c1995$cpi-c1985$tax/c1985$cpi,
                   dInc=log(c1995$rincome/c1985$rincome))
res1 <- tsls(dQ~dP+dInc, ~dInc+dTs, vcov="MDS", data=data)
res2 <- tsls(dQ~dP+dInc, ~dInc+dT, vcov="MDS", data=data)
res3 <- tsls(dQ~dP+dInc, ~dInc+dTs+dT, vcov="MDS", data=data)

## -----------------------------------------------------------------------------
res4 <- gmm4(dQ~dP+dInc, ~dInc+dTs+dT, vcov="MDS", data=data)
specTest(res4)

## ----echo=FALSE, results='asis'-----------------------------------------------
texreg(list(res1,res2,res3), fontsize="footnotesize", label="tab1", 
       caption="Table 12.1 of Stock and Watson textbook", df.adj=TRUE, sandwich=TRUE)

## -----------------------------------------------------------------------------
res4 <- gmm4(dQ~dP+dInc, ~dInc+dTs+dT, vcov="MDS", data=data)

## -----------------------------------------------------------------------------
data(HealthRWM)
dat88 <- subset(HealthRWM, year==1988 & hhninc>0)
dat88$hhninc <- dat88$hhninc/10000

## -----------------------------------------------------------------------------
thet0 <- c(b0=log(mean(dat88$hhninc)),b1=0,b2=0,b3=0)
g <- hhninc~exp(b0+b1*age+b2*educ+b3*female)
res0 <- nls(g, dat88, start=thet0, control=list(maxiter=100))
summary(res0)$coef

## -----------------------------------------------------------------------------
h1 <- ~age+educ+female
model1 <- momentModel(g, h1, thet0, vcov="MDS", data=dat88)
res1 <- gmmFit(model1, control=list(reltol=1e-10, abstol=1e-10))

## -----------------------------------------------------------------------------
h2 <- ~age+educ+female+hsat+married
model2 <- momentModel(g, h2, thet0, vcov="MDS", data=dat88)
res2 <- gmmFit(model2, type="onestep")

## -----------------------------------------------------------------------------
res3 <- gmmFit(model2)

## ----echo=FALSE, results='asis'-----------------------------------------------
texreg(list(res1, res2, res3), caption="Attempt to reproduce Table 13.2 from Greene (2012)",
       label="greene1", fontsize='footnotesize', digits=5,
       includeJTest=FALSE, includeFTest=FALSE)

## -----------------------------------------------------------------------------
data(ConsumptionG)
Y <- ConsumptionG$REALDPI
C <- ConsumptionG$REALCONS
n <- nrow(ConsumptionG)
Y1 <- Y[-n]; Y <- Y[-1]
C1 <- C[-n]; C <- C[-1]
dat <- data.frame(Y=Y,Y1=Y1,C=C,C1=C1)
model <- momentModel(C~Y, ~Y1+C1, data=dat, vcov="iid")

## -----------------------------------------------------------------------------
res1 <- tsls(model)
res2 <- lm(C~Y)

## -----------------------------------------------------------------------------
DWH(res1)

## -----------------------------------------------------------------------------
DWH(res1, res2, df.adj=TRUE)

## -----------------------------------------------------------------------------
X <- model.matrix(model)
Xhat <- qr.fitted(res1@wObj@w, X)
s2 <- sum(residuals(res2)^2)/(res2$df.residual)
v1 <-  solve(crossprod(Xhat))*s2
v2 <- solve(crossprod(X))*s2
DWH(res1, res2, v1=v1, v2=v2)

## -----------------------------------------------------------------------------
data(simData)
g <- list(Supply=y1~x1+z2, Demand1=y2~x1+x2+x3, Demand2=y3~x3+x4+z1)
h <- list(~z1+z2+z3, ~x3+z1+z2+z3+z4, ~x3+x4+z1+z2+z3)
smod1 <- sysMomentModel(g, h, vcov="iid", data=simData)
smod1

## -----------------------------------------------------------------------------
smod2 <- sysMomentModel(g, ~x2+x4+z1+z2+z3+z4, vcov="iid", data=simData)
smod2

## -----------------------------------------------------------------------------
smod3 <- sysMomentModel(g, vcov="iid", data=simData)
smod3

## -----------------------------------------------------------------------------
dat <- list(y=matrix(rnorm(150),50,3),
            x=rnorm(50), z1=rnorm(50),
            z2=rnorm(50))
mod <- momentModel(y~x, ~z1+z2, data=dat, vcov="iid")
mod

## -----------------------------------------------------------------------------
mod <- momentModel(y~x, ~x, vcov="iid", data=dat)

## -----------------------------------------------------------------------------
setCoef(smod1, 1:11)

## -----------------------------------------------------------------------------
smod1[1:2, list(1:3,1:4)]

## -----------------------------------------------------------------------------
gmmFit(smod1[1])

## -----------------------------------------------------------------------------
mm <- model.matrix(smod1)
mm <- lapply(1:3, function(i) model.matrix(smod1[i]))

## -----------------------------------------------------------------------------
theta <- list(1:3, 1:4, 1:4)
gt <- evalMoment(smod1, theta)

## -----------------------------------------------------------------------------
Sigma <- crossprod(residuals(smod1, theta))/smod1@n

## -----------------------------------------------------------------------------
eq1 <- momentModel(g[[1]], h[[1]], data=simData, vcov="iid")
eq2 <- momentModel(g[[2]], h[[2]], data=simData, vcov="iid")
eq3 <- momentModel(g[[3]], h[[3]], data=simData, vcov="iid")
smod <- merge(eq1,eq2,eq3)
smod

## -----------------------------------------------------------------------------
eq1 <- momentModel(y~x1, ~x1+z4, data=simData, vcov="iid")
merge(smod1, eq1)

## -----------------------------------------------------------------------------
R1 <- list(c("x1=-12*z2"), character(), c("x3=0.8", "z1=0.3"))
rsmod1 <- restModel(smod1, R1)
rsmod1

## -----------------------------------------------------------------------------
R2<- c("Supply.x1=1", "Demand1.x3=Demand2.x3")
rsmod1.ce <- restModel(smod1, R2)
rsmod1.ce

## -----------------------------------------------------------------------------
e <- residuals(rsmod1, theta=list(1:2, 1:4, 1:2))
dim(e)

## -----------------------------------------------------------------------------
(b <- coef(rsmod1, theta=list(1:2, 1:4, 1:2)))
e <- residuals(as(rsmod1, "slinearModel"), b)

## -----------------------------------------------------------------------------
evalDMoment(rsmod1, theta=list(1:2,1:4,1:2))[[1]]

## -----------------------------------------------------------------------------
e <- residuals(rsmod1.ce, theta=list(1:9))
e[1:3,]

## -----------------------------------------------------------------------------
(b <- coef(rsmod1.ce, theta = list(1:9)))
e <- residuals(as(rsmod1.ce, "slinearModel"), b)

## -----------------------------------------------------------------------------
G <- evalDMoment(rsmod1.ce, list(1:9))
names(G)
dim(G[[1]])

## -----------------------------------------------------------------------------
rsmod1[1]

## -----------------------------------------------------------------------------
nsmod <- as(smod1, "snonlinearModel")
nsmod

## -----------------------------------------------------------------------------
nsmod@parNames

## -----------------------------------------------------------------------------
setCoef(nsmod, 1:11)

## -----------------------------------------------------------------------------
R1 <- c("theta1=-12*theta2","theta9=0.8", "theta11=0.3")
R2<- c("theta1=1", "theta6=theta10")
(rnsmod1 <- restModel(nsmod, R1))
(rnsmod2 <- restModel(nsmod, R2))

## -----------------------------------------------------------------------------
wObj1 <- evalWeights(smod1, w="ident")
wObj1

## -----------------------------------------------------------------------------
wObj1@sameMom

## -----------------------------------------------------------------------------
wObj1@type

## -----------------------------------------------------------------------------
wObj1@eqnNames
wObj1@momNames

## -----------------------------------------------------------------------------
wObj2 <- evalWeights(smod1, w=diag(16))

## -----------------------------------------------------------------------------
smod1 <- sysMomentModel(g,h,vcov="MDS", data=simData)
wObj <- evalWeights(smod1, theta=list(1:3,1:4,1:4))
is(wObj@w)
wObj@Sigma

## -----------------------------------------------------------------------------
gt <- evalMoment(smod1, theta=list(1:3, 1:4, 1:4)) ## this is a list
gbar <- colMeans(do.call(cbind, gt))
obj <- smod1@n*quadra(wObj, gbar)
obj

## -----------------------------------------------------------------------------
evalGmmObj(smod1, theta=list(1:3,1:4,1:4), wObj=wObj)

## -----------------------------------------------------------------------------
smod1 <- sysMomentModel(g,h,vcov="MDS", data=simData)
wObj1 <- evalWeights(smod1, w="ident")
theta0 <- solveGmm(smod1, wObj1)$theta
wObj2 <- evalWeights(smod1, theta=theta0)
solveGmm(smod1, wObj2)

## -----------------------------------------------------------------------------
R1 <- list(c("x1=-12*z2"), character(), c("x3=0.8", "z1=0.3"))
rsmod1 <- restModel(smod1, R1)
wObj1 <- evalWeights(rsmod1, w="ident")
theta0 <- solveGmm(rsmod1, wObj1)$theta
wObj2 <- evalWeights(rsmod1, theta=theta0)
theta1 <- solveGmm(rsmod1, wObj2)$theta
theta1

## -----------------------------------------------------------------------------
coef(rsmod1, theta1)

## -----------------------------------------------------------------------------
R2<- c("Supply.x1=1", "Demand1.x3=Demand2.x3")
rsmod1<- restModel(smod1, R2)
wObj1 <- evalWeights(rsmod1, w="ident")
theta0 <- solveGmm(rsmod1, wObj1)$theta
wObj2 <- evalWeights(rsmod1, theta=theta0)
theta1 <- solveGmm(rsmod1, wObj2)$theta
theta1

## -----------------------------------------------------------------------------
coef(rsmod1, theta1)

## -----------------------------------------------------------------------------
### Without cross-equation restrictions
wObj1 <- evalWeights(rnsmod1, w="ident")
theta0 <- solveGmm(rnsmod1, wObj1, theta0=rep(0, 8))$theta
wObj2 <- evalWeights(rnsmod1, theta=theta0)
theta1 <- solveGmm(rnsmod1, wObj2, theta0=theta0)$theta
### Verify that the restrictions are correctly imposed:
printRestrict(rnsmod1)
coef(rnsmod1, theta1)
### With cross-equation restrictions
wObj1 <- evalWeights(rnsmod2, w="ident")
theta0 <- solveGmm(rnsmod2, wObj1, theta0=rep(0, 9))$theta
wObj2 <- evalWeights(rnsmod2, theta=theta0)
theta1 <- solveGmm(rnsmod2, wObj2, theta0=theta0)$theta
### Verify that the restrictions are correctly imposed:
printRestrict(rnsmod2)
coef(rnsmod2, theta1)

## -----------------------------------------------------------------------------
smod1 <- sysMomentModel(g,h,vcov="MDS", data=simData)
gmmFit(smod1, type="twostep")

## -----------------------------------------------------------------------------
smod1 <- sysMomentModel(g,h,vcov="iid", data=simData)
gmmFit(smod1, type="twostep")

## -----------------------------------------------------------------------------
smod1 <- sysMomentModel(g,~z1+z2+z3+z4+z5,vcov="iid", data=simData)
gmmFit(smod1, type="twostep", initW="tsls")

## -----------------------------------------------------------------------------
smod1 <- sysMomentModel(g, vcov="iid", data=simData)
gmmFit(smod1, type="twostep", initW="tsls")

## -----------------------------------------------------------------------------
smod1 <- sysMomentModel(g,h,vcov="MDS", data=simData)
res <- gmmFit(smod1, type="twostep", initW="EbyE")

## -----------------------------------------------------------------------------
gmmFit(smod1,  EbyE=TRUE) ## type is 'twostep' by default

## -----------------------------------------------------------------------------
res <- gmmFit(smod1,  EbyE=TRUE, weights="ident")

## -----------------------------------------------------------------------------
R1 <- list(c("x1=-12*z2"), character(), c("x3=0.8", "z1=0.3"))
rsmod1 <- restModel(smod1, R1)
gmmFit(rsmod1)@theta
R2<- c("Supply.x1=1", "Demand1.x3=Demand2.x3")
rsmod1<- restModel(smod1, R2)
gmmFit(rsmod1)@theta

## -----------------------------------------------------------------------------
theta0 <- setCoef(rnsmod1, rep(0,8))
gmmFit(rnsmod1, theta0=theta0)@theta
theta0 <- setCoef(rnsmod2, rep(0,9))
gmmFit(rnsmod2, theta0=theta0)@theta

## -----------------------------------------------------------------------------
smod1 <- sysMomentModel(g,h,vcov="MDS", data=simData)
res <- tsls(smod1)
res

## -----------------------------------------------------------------------------
smod2 <- sysMomentModel(g,~z1+z2+z3+z4+z5,vcov="MDS", data=simData)
res <- ThreeSLS(smod2)

## -----------------------------------------------------------------------------
smod2 <- sysMomentModel(g,,vcov="MDS", data=simData)
res <- ThreeSLS(smod2)

## -----------------------------------------------------------------------------
smod2 <- sysMomentModel(g,~z1+z2+z3+z4+z5,vcov="iid", data=simData)
gmmFit(smod2, initW="tsls")@theta
ThreeSLS(smod2)@theta

## -----------------------------------------------------------------------------
smod1 <- sysMomentModel(g, h, vcov="iid", data=simData)
res <- gmmFit(smod1)
specTest(res)

## -----------------------------------------------------------------------------
summary(res)

## -----------------------------------------------------------------------------
smod1 <- sysMomentModel(g,h,vcov="iid", data=simData)
R1 <- list(c("x1=-12*z2"), character(), c("x3=0.8", "z1=0.3"))
rsmod1 <- restModel(smod1, R1)
summary(gmmFit(rsmod1))@coef
R2<- c("Supply.x1=1", "Demand1.x3=Demand2.x3")
rsmod1<- restModel(smod1, R2)
summary(gmmFit(rsmod1))@coef

## -----------------------------------------------------------------------------
smod1 <- sysMomentModel(g, h, vcov="MDS", data=simData)
res.u <- gmmFit(smod1)
R1 <- list(c("x1=-12*z2"), character(), c("x3=0.8", "z1=0.3"))
rsmod1 <- restModel(smod1, R1)
res.r <- gmmFit(rsmod1)

## -----------------------------------------------------------------------------
hypothesisTest(res.u, res.r, type="Wald")

## -----------------------------------------------------------------------------
R2<- c("Supply.x1=1", "Demand1.x3=Demand2.x3")
rsmod1<- restModel(smod1, R2)
res2.r <- gmmFit(rsmod1)
hypothesisTest(res.u, res2.r, type="LR")

## -----------------------------------------------------------------------------
R1 <- c("theta1=-12*theta2","theta9=0.8", "theta11=0.3")
R2<- c("theta1=1", "theta6=theta10")
rnsmod1 <- restModel(nsmod, R1)
rnsmod2 <- restModel(nsmod, R2)
theta0 <- setCoef(nsmod, rep(0,11))
fit <- gmmFit(nsmod, theta0=theta0)
theta0 <- setCoef(rnsmod1, rep(0,8))
rfit1 <- gmmFit(rnsmod1, theta0=theta0)
theta0 <- setCoef(rnsmod2, rep(0,9))
rfit2 <- gmmFit(rnsmod2, theta0=theta0)

## -----------------------------------------------------------------------------
hypothesisTest(object.u=fit, R=R1)
hypothesisTest(object.u=fit, object.r=rfit1, type="LR")
hypothesisTest(object.u=fit, object.r=rfit1, type="LM")
hypothesisTest(object.u=fit, R=R2)
hypothesisTest(object.u=fit, object.r=rfit2, type="LR")
hypothesisTest(object.u=fit, object.r=rfit2, type="LM")

## -----------------------------------------------------------------------------
res <- gmm4(g, h, type="twostep", vcov="MDS", data=simData)
res

## -----------------------------------------------------------------------------
res <- gmm4(g, h, type="twostep", vcov="MDS", EbyE=TRUE, data=simData)
res

## -----------------------------------------------------------------------------
res <- gmm4(g, ~z1+z2+z3+z4+z5, type="twostep", vcov="iid", initW="tsls", data=simData) #3SLS
res <- gmm4(g, NULL, type="twostep", vcov="iid", initW="tsls", data=simData) #SUR

## -----------------------------------------------------------------------------
R1 <- list(c("x1=-12*z2"), character(), c("x3=0.8", "z1=0.3"))
res <- gmm4(g, h, data=simData, cstLHS=R1) #two-step by default
res

## -----------------------------------------------------------------------------
h <- list(~z1+z2+z3, ~x3+z1+z2+z3+z4, ~x3+x4+z1+z2+z3)
nlg <- list(Supply=y1~theta0+theta1*x1+theta2*z2,
            Demand1=y2~alpha0+alpha1*x1+alpha2*x2+alpha3*x3,
            Demand2=y3~beta0+beta1*x3+beta2*x4+beta3*z1)
theta0 <- list(c(theta0=0,theta1=0,theta2=0),
               c(alpha0=0,alpha1=0,alpha2=0, alpha3=0),
               c(beta0=0,beta1=0,beta2=0,beta3=0))
fit <- gmm4(nlg, h, theta0,data=simData)
## the restricted estimation (:
R2<- c("theta1=1", "alpha1=beta2")
fit2 <- gmm4(nlg, h, theta0,data=simData, cstLHS=R2)

## -----------------------------------------------------------------------------
data(ManufactCost)
price <- c("Pk","Pl","Pe")
ManufactCost[,price] <- log(ManufactCost[,price]/ManufactCost$Pm)

## -----------------------------------------------------------------------------
g <- list(Sk=K~Pk+Pl+Pe,
          Sl=L~Pk+Pl+Pe,
          Se=E~Pk+Pl+Pe)
mod <- sysMomentModel(g, NULL, data=ManufactCost, vcov="iid")

## -----------------------------------------------------------------------------
R <- c("Sk.Pl=Sl.Pk", "Sk.Pe=Se.Pk", "Sl.Pe=Se.Pl")
rmod <- restModel(mod, R=R)

## -----------------------------------------------------------------------------
res <- gmmFit(rmod)
summary(res)@coef

## -----------------------------------------------------------------------------
res.u <- gmmFit(mod)
hypothesisTest(res.u, res)

## -----------------------------------------------------------------------------
data(Klein)
Klein1 <- Klein[-22,]
Klein <- Klein[-1,]
dimnames(Klein1) <- list(rownames(Klein), paste(colnames(Klein),"1",sep=""))
Klein <- cbind(Klein, Klein1)
Klein$A <- (Klein$YEAR-1931)

## -----------------------------------------------------------------------------
g <- list(C=C~P+P1+I(WP+WG),
          I=I~P+P1+K1,
          Wp=WP~X+X1+A)
h <- ~G+T+WG+A+K1+P1+X1          
res <- ThreeSLS(g, h, vcov="iid", data=Klein)
summary(res, breadOnly=TRUE)@coef

## ----extract, eval=FALSE------------------------------------------------------
#  library(texreg)
#  setMethod("extract", "gmmfit",
#            function(model, includeJTest=TRUE, includeFTest=TRUE, ...)
#                {
#                    s <- summary(model, ...)
#                    spec <- modelDims(model@model)
#                    coefs <- s@coef
#                    names <- rownames(coefs)
#                    coef <- coefs[, 1]
#                    se <- coefs[, 2]
#                    pval <- coefs[, 4]
#                    n <- model@model@n
#                    gof <- numeric()
#                    gof.names <- character()
#                    gof.decimal <- logical()
#                    if (includeJTest) {
#                        if (spec$k == spec$q)
#                            {
#                                obj.fcn <- NA
#                                obj.pv <- NA
#                            } else {
#                                obj.fcn <- s@specTest@test[1]
#                                obj.pv <- s@specTest@test[3]
#                            }
#                        gof <- c(gof, obj.fcn, obj.pv)
#                        gof.names <- c(gof.names, "J-test Statistics", "J-test p-value")
#                        gof.decimal <- c(gof.decimal, TRUE, TRUE)
#                    }
#                    if (includeFTest) {
#                        str <- s@strength$strength
#                        if (is.null(str))
#                            {
#                                gof <- c(gof, NA)
#                                gof.names <- c(gof.names, "First Stage F-stats")
#                                gof.decimal <- c(gof.decimal, TRUE)
#                            } else {
#                                for (i in 1:nrow(str))
#                                    {
#                                        gof <- c(gof, str[i,1])
#                                        gofn <- paste("First Stage F-stats(",
#                                                      rownames(str)[i], ")", sep="")
#                                        gof.names <- c(gof.names, gofn)
#                                        gof.decimal <- c(gof.decimal, TRUE)
#                                    }
#                            }
#                    }
#                    tr <- createTexreg(coef.names = names, coef = coef, se = se,
#                                       pvalues = pval, gof.names = gof.names, gof = gof,
#                                       gof.decimal = gof.decimal)
#                    return(tr)
#                })

