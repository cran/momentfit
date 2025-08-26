## ----echo=FALSE---------------------------------------------------------------
library(knitr)
opts_chunk$set(size='footnotesize')

## -----------------------------------------------------------------------------
library(ivmodel)
data(card.data)
## from the ivmodel examples
Z <- card.data[,c("nearc4","nearc2")]
Y <- card.data[,"lwage"]
D <- card.data[,"educ"]
Xname <- c("exper", "expersq", "black", "south", "smsa", "reg661",
           "reg662", "reg663", "reg664", "reg665", "reg666", "reg667",
           "reg668", "smsa66")
X <- card.data[,Xname]
mod <- ivmodel(Y=Y,D=D,Z=Z,X=X)

## -----------------------------------------------------------------------------
c(LIML=mod$LIML$k, Fuller=mod$Fuller$k)

## -----------------------------------------------------------------------------
library(momentfit)
g <- reformulate(c("educ", Xname), "lwage")
h <- reformulate(c(c("nearc4","nearc2"), Xname))
mod2 <- momentModel(g, h, data=card.data, vcov="iid")

## -----------------------------------------------------------------------------
getK(mod2)

## -----------------------------------------------------------------------------
g2 <- reformulate(c("educ", "educ:exper", Xname), "lwage")
h2 <- reformulate(c(c("nearc4","nearc2", "nearc2:exper", "nearc4:exper"), Xname))
mod3 <- momentModel(g2, h2, data=card.data)
getK(mod3)

## -----------------------------------------------------------------------------
h3 <- reformulate(c(c("nearc4"), Xname))
mod4 <- momentModel(g, h3, data=card.data)
getK(mod4)

## -----------------------------------------------------------------------------
lse(mod2)

## -----------------------------------------------------------------------------
(liml <- kclassfit(mod2))
(fuller <- kclassfit(mod2, type="Fuller"))
(btsls <- kclassfit(mod2, type="BTSLS"))

## -----------------------------------------------------------------------------
print(mod$LIML$point.est,digits=10)
print(coef(liml)[2], digits=10)

## -----------------------------------------------------------------------------
print(mod$Fuller$point.est,digits=10)
print(coef(fuller)[2], digits=10)

## -----------------------------------------------------------------------------
resK <- getK(mod2, 1, TRUE)
(liml <- kclassfit(mod2, resK))
(fuller <- kclassfit(mod2, resK, type="Fuller"))

## -----------------------------------------------------------------------------
(s <- summary(liml))

## -----------------------------------------------------------------------------
specTest(liml)

## -----------------------------------------------------------------------------
s@coef["educ",]
mod$LIML$std.err

## -----------------------------------------------------------------------------
spec <- modelDims(mod2)
u <- residuals(liml)
sig <- sum(u^2)/(spec$n-spec$k)
W <- model.matrix(liml@model, "instruments")
myX <- model.matrix(liml@model)
sqrt(diag(sig*solve(t(W)%*%myX)))[2]

## ----eval=FALSE---------------------------------------------------------------
# mod <- ivmodel(Y=Y,D=D,Z=Z,X=X,heteroSE=TRUE)
# mod2 <- momentModel(g, h, data=card.data, vcov="MDS")
# liml <- kclassfit(mod2, resK)
# summary(liml)@coef["educ",]
# c(mod$LIML$point.est, mod$LIML$std.err)

## -----------------------------------------------------------------------------
(CD2 <- CDtest(mod2, print=FALSE))
momentStrength(mod2)

## -----------------------------------------------------------------------------
CDtest(mod2)

## -----------------------------------------------------------------------------
g <- reformulate(c("educ", Xname), "lwage")
h <- reformulate(c(c("nearc4","nearc2","IQ","KWW"), Xname))
mod5 <- momentModel(g, h, data=card.data, vcov="iid")
CDtest(mod5)

## ----message=FALSE------------------------------------------------------------
data(simData)
## Step 1
m <- tsls(y1~y2+y3+x1+x2, ~z1+z2+z3+z4+z5+x1+x2, data=simData)
e <- residuals(m)
## Step 2
fit <- lm(e~z1+z2+z3+z4+z5+x1+x2, simData)
fitr <- lm(e~x1+x2, simData)
F <- anova(fit, fitr)$F[2]
## Step 4
(sw1 <- F*5/(5-2))

## -----------------------------------------------------------------------------
smod <- momentModel(y~y1+y2+y3+x1+x2, ~z1+z2+z3+z4+z5+x1+x2, data=simData)
SWtest(smod,1,FALSE)

## -----------------------------------------------------------------------------
SWtest(smod)

## -----------------------------------------------------------------------------
MOPtest(mod2)

## -----------------------------------------------------------------------------
CDtest(mod2, FALSE)

## -----------------------------------------------------------------------------
mod2 <- update(mod2, vcov="MDS")
MOPtest(mod2)

## -----------------------------------------------------------------------------
## get Z
spec <- modelDims(mod2)
Z2 <- model.matrix(mod2, "excludedExo")
X1 <- model.matrix(mod2, "includedExo")
X2 <- model.matrix(mod2, "includedEndo")
y <- modelResponse(mod2)

## -----------------------------------------------------------------------------
fitX1  <- lm.fit(X1, Z2)
Z2 <- fitX1$residuals
X2 <- qr.resid(fitX1$qr, X2)
y <- qr.resid(fitX1$qr, y)
Z2 <- qr.Q(qr(Z2))*sqrt(nrow(Z2))

## -----------------------------------------------------------------------------
colnames(Z2) = paste("Z", 1:ncol(Z2), sep="")
dat <- as.data.frame(cbind(X2,Z2))
g <- reformulate(colnames(Z2), colnames(X2), FALSE)
h <- reformulate(colnames(Z2), NULL, FALSE)
m2 <- momentModel(g,h,data=dat,vcov="MDS")

## -----------------------------------------------------------------------------
b <- crossprod(Z2,X2)/nrow(Z2)
w2 <- vcov(m2, b)

## -----------------------------------------------------------------------------
dat <- as.data.frame(cbind(y=y,X2,Z2))
g <- list(reformulate(colnames(Z2), "y", FALSE),  g)
h <- list(h)
m <- sysMomentModel(g,h,data=dat,vcov="MDS")
b <- list(c(crossprod(Z2,y)/nrow(Z2)),
          c(crossprod(Z2,X2)/nrow(Z2)))
w <- vcov(m, b)
w

## -----------------------------------------------------------------------------
w2

## -----------------------------------------------------------------------------
res <- momentfit:::getMOPw(mod2)
res$w2

## -----------------------------------------------------------------------------
momentfit:::getMOPx(w=res, tau=0.10, type="LIML")

## -----------------------------------------------------------------------------
MOPtest(mod5, simplified=FALSE, estMethod="LIML")

## -----------------------------------------------------------------------------
MOPtest(mod5,  simplified=FALSE, estMethod="TSLS")

## -----------------------------------------------------------------------------
mod4 <- update(mod4, vcov="MDS")
MOPtest(mod4, estMethod="TSLS", print=FALSE)["Feff"]

## -----------------------------------------------------------------------------
momentStrength(update(mod4, vcovOptions=list(type="HC0")))$strength

## -----------------------------------------------------------------------------
LewMertest(mod3, simplified=FALSE)

## -----------------------------------------------------------------------------
LewMertest(mod3, simplified=FALSE, print=FALSE, npoints=1)$crit
LewMertest(mod3, simplified=FALSE, print=FALSE, npoints=20)$crit

## -----------------------------------------------------------------------------
getIVDat <- function(n, R2, k, rho, b0=0)
{
    eta <- sqrt(R2/(k*(1-R2)))
    Z <- sapply(1:k, function(i) rnorm(n))
    sigma <- chol(matrix(c(1,rho,rho,1),2,2))
    err <- cbind(rnorm(n), rnorm(n))%*%sigma
    y2 <- rowSums(Z)*eta+err[,2]
    y1 <- b0*y2 + err[,1]
    dat <- data.frame(y1=y1, y2=y2, u=err[,1], e=err[,2])
    for (i in 1:k) dat[[paste("Z",i,sep="")]] <- Z[,i]
    dat
}

## -----------------------------------------------------------------------------
library(momentfit)
set.seed(112233)
k <- 10
rho <- .3
R2 <- .001
g <- y1~y2
n <- 500
h <- reformulate(paste("Z", 1:k, sep=""))
dat <- getIVDat(n, R2, k, rho)
m <- momentModel(g, h, data=dat, vcov="MDS")

## ----echo=FALSE---------------------------------------------------------------
## n: sample size
## N: number of endogenous
## K: number of Instruments
## Nx: number of exogenous regressors

genDat <- function(n=200, N=2, K=5, Nx=2)
{
    beta <- rep(1,N) 
    DL <- 0.08
    X <- qr.R(qr(matrix(rnorm(K^2),K,K)))
    L0 <- t(X[,1:N])
    Pi <- sqrt(K)*DL^0.5*L0
    Pi <- t(Pi)
    By <- rnorm(Nx+1)
    BY <- matrix(rnorm((Nx+1)*N), Nx+1, N)
    u <- rnorm(n)
    v <- matrix(rnorm(n*N), n, N)
    X <- cbind(matrix(rnorm(n*Nx), n, Nx), 1)
    Z <- matrix(rnorm(n*K), n, K)
    Y <- Z%*%Pi+X%*%BY+v
    y <- c(Y%*%beta+X%*%By+u)
    X <- X[,-ncol(X),drop=FALSE]    
    colnames(Y) <- paste("Y", 1:ncol(Y), sep="")
    colnames(X) <- paste("X", 1:ncol(X), sep="")
    colnames(Z) <- paste("Z", 1:ncol(Z), sep="")
    dat <- as.data.frame(cbind(Y,X,Z))
    dat$y <- y
    dat
}
dat <- genDat()
m <- momentModel(y~Y1+X1+X2, ~Z1+Z2+Z3+Z4+Z5+X1+X2, data=dat,
                 vcov="MDS")
LewMertest(m, simplified=FALSE, npoints=100)
MOPtest(m, estMethod="TSLS", simplified=FALSE)

