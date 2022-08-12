#####  All S4 classes of the package are defined here
######################################################


## Union Classes

setClassUnion("matrixORcharacter", c("matrix", "character"))
setClassUnion("matrixORnumeric", c("matrix", "numeric"))
setClassUnion("numericORcharacter", c("numeric", "character"))
setClassUnion("numericORNULL", c("numeric", "NULL"))
setClassUnion("numericORlogical", c("numeric", "logical"))
setClassUnion("numericORmatrixORNULL", c("matrix", "numeric", "NULL"))
setClassUnion("expressionORNULL", c("expression", "NULL"))
setClassUnion("functionORNULL", c("function", "NULL"))
setClassUnion("callORNULL", c("call", "NULL"))


## smooth spec class

setOldClass("tskernel")
setClass("sSpec", slots=list(k="numeric", kernel="character", bw="numeric",w="tskernel",
                             bwMet="character"),
         prototype=list(w=kernel(1), bw=1, k=c(1,1), kernel="none", bwMet="none"))

## moment based models
setClass("linearModel", slots = list(modelF="data.frame", instF="data.frame",
                                     vcov="character",n="integer", q="integer", k="integer",
                                     parNames="character", momNames="character",
                                     vcovOptions="list", centeredVcov="logical",
                                     varNames="character", isEndo="logical",
                                     omit='integer', survOptions="list",
                                     sSpec="sSpec", smooth="logical"))
setClass("nonlinearModel", slots = list(modelF="data.frame", instF="data.frame",
                                        vcov="character",theta0="numeric",
                                        n="integer", q="integer",k="integer",
                                        parNames="character", momNames="character",
                                        fRHS="expression", fLHS="expressionORNULL",
                                        vcovOptions="list",
                                        centeredVcov="logical", varNames="character",
                                        isEndo="logical",omit='integer', survOptions="list",
                                        sSpec="sSpec", smooth="logical"))
setClass("functionModel", slots = list(X="ANY", fct="function",dfct="functionORNULL",
                                       vcov="character",theta0="numeric",
                                       n="integer", q="integer",k="integer",
                                       parNames="character", momNames="character",
                                       vcovOptions="list",
                                       centeredVcov="logical", varNames="character",
                                       isEndo="logical",omit='integer', survOptions="list",
                                       sSpec="sSpec", smooth="logical"))
setClass("formulaModel", slots = list(modelF="data.frame", 
                                        vcov="character",theta0="numeric",
                                        n="integer", q="integer",k="integer",
                                        parNames="character", momNames="character",
                                        fRHS="list", fLHS="list",
                                        vcovOptions="list",
                                        centeredVcov="logical", varNames="character",
                                        isEndo="logical", isMDE="logical",omit='integer',
                                        survOptions="list",sSpec="sSpec", smooth="logical"))
setClassUnion("regModel", c("linearModel", "nonlinearModel"))
setClassUnion("allNLModel", c("nonlinearModel", "functionModel", "formulaModel"))
setClassUnion("momentModel", c("linearModel", "nonlinearModel", "functionModel", "formulaModel"))

## Restricted Models

setClass("rlinearModel", slots = list(cstLHS="matrix", cstRHS="numeric", cstSpec="list"),
         contains="linearModel")

setClass("rnonlinearModel", slots = list(R="list", cstSpec="list"),
         contains="nonlinearModel")

setClass("rfunctionModel", slots = list(R="list", cstSpec="list"),
         contains="functionModel")

setClass("rformulaModel", slots = list(R="list", cstSpec="list"),
         contains="formulaModel")

setClassUnion("rmomentModel", c("rlinearModel", "rnonlinearModel", "rfunctionModel",
                                "rformulaModel"))

### System models

setClass("slinearModel", slots = list(modelT="list", instT="list",data="data.frame",
                                      vcov="character",n="integer", q="integer",
                                      k="integer", parNames="list",
                                      momNames="list", eqnNames="character",
                                      vcovOptions="list",
                                      centeredVcov="logical", sameMom="logical",
                                      SUR="logical", varNames="list", isEndo="list",
                                      omit='integer', survOptions="list",
                                      sSpec="sSpec", smooth="logical"))

setClass("snonlinearModel", slots = list(data="data.frame", instT="list",
                                         vcov="character",theta0="list",
                                         n="integer", q="integer",k="integer",
                                         parNames="list", momNames="list",
                                         fRHS="list", fLHS="list", eqnNames="character",
                                         vcovOptions="list",
                                         centeredVcov="logical", sameMom="logical",
                                         SUR="logical",
                                         varNames="list", isEndo="list",
                                         omit='integer', survOptions="list",
                                         sSpec="sSpec", smooth="logical"))

setClass("sfunctionModel", slots = list(X="ANY", fct="list", dfct="list",
                                        vcov="character",theta0="list",
                                        n="integer", q="integer",k="integer",
                                        parNames="list", momNames="list",
                                        eqnNames="character", vcovOptions="list",
                                        centeredVcov="logical", 
                                        varNames="list",
                                        sameMom="logical", SUR="logical",
                                        omit='integer', survOptions="list",
                                        sSpec="sSpec", smooth="logical"))

setClassUnion("sysModel", c("slinearModel", "snonlinearModel",
                            "sfunctionModel"))



## Restricted System models

setClass("rslinearModel", slots = list(cstLHS="matrix", cstRHS="numeric", cstSpec="list"),
         contains="slinearModel")

setClass("rsnonlinearModel", slots = list(R="list", cstSpec="list"),
         contains="snonlinearModel")

setClassUnion("rsysModel", c("rslinearModel", "rsnonlinearModel"))



## momentWeights

setClass("momentWeights", representation(w="ANY", type="character", wSpec="list"))

### sysMomentWeights

setClass("sysMomentWeights", representation(w="ANY", type="character", wSpec="list",
                                         Sigma="ANY", momNames="list", eqnNames="character",
                                         sameMom="logical"))


## specTest

setClass("specTest", representation(test = "matrix", testname="character"))

## gmmfit

setClass("gmmfit", representation(theta = "numeric", convergence = "numericORNULL",
                                  convIter="numericORNULL",call="callORNULL",
                                  type="character", wObj="momentWeights",niter="integer",
                                  efficientGmm="logical", model="momentModel"))

## summaryGmm


setClass("summaryGmm", representation(coef="matrix", specTest = "specTest",
                                      strength="list", model="momentModel",sandwich="logical",
                                      type="character", convergence = "numericORNULL",
                                      convIter="numericORNULL", wSpec="list",niter="integer",
                                      df.adj="logical", breadOnly="logical"))


## hypothesisTest

setClass("hypothesisTest", representation(test="numeric", hypothesis="character",
                                          dist="character", df="integer", pvalue="numeric",
                                          type="character"))


## summarySysGmm

setClass("summarySysGmm",
         representation(coef="list", specTest = "specTest",
                        strength="list", model="sysModel",sandwich="logical",
                        type="character", convergence = "numericORNULL",
                        convIter="numericORNULL", wSpec="list",niter="integer",
                        df.adj="logical", breadOnly="logical"))


## "tsls"

setClass("tsls", contains="gmmfit")

## 

## confint

setClass("confint", representation(interval = "matrix", type="character",
                                   level="numeric", theta="numeric"))


setClass("mconfint", 
         representation(areaPoints="matrix", type="character", level="numeric",
                        theta="numeric"))


### system GMM fit

setClass("sgmmfit", representation(theta = "list", convergence = "numericORNULL",
                                   convIter="numericORNULL",call="callORNULL",
                                   type="character", wObj="sysMomentWeights",niter="integer",
                                   efficientGmm="logical", model="sysModel"))

setClass("stsls", contains="sgmmfit")

## gelfit

setClass("gelfit", representation(theta = "numeric", convergence = "numeric",
                                  lambda = "numeric", lconvergence = "numeric",
                                  call="callORNULL", gelType="list", vcov="list",
                                  model="momentModel", restrictedLam="integer"))

setClass("summaryGel", representation(coef="matrix", specTest = "specTest",
                                      model="momentModel", lambda="matrix",
                                      convergence="numeric",lconvergence="numeric",
                                      impProb="list", gelType="list",
                                      restrictedLam="integer"))


## class converted

setAs("linearModel", "nonlinearModel",
      function(from) {
          spec <- modelDims(from)
          X <- model.matrix(from)
          theta0 <- rep(1,ncol(X))
          names(theta0) <- paste("theta", 1:ncol(X), sep="")         
          colnames(X) <- paste("X", 1:ncol(X), sep="")
          rhs <- paste(names(theta0), "*", colnames(X), sep="")
          rhs <- paste(rhs, collapse="+", sep="")
          rhs <- parse(text=rhs)
          X <- data.frame(Y=modelResponse(from), X)
          lhs <- expression(Y)
          new("nonlinearModel", modelF=X, instF=from@instF, vcov=from@vcov,
              theta0=theta0, n=spec$n, q=spec$q, k=spec$k, parNames=names(theta0),
              momNames=spec$momNames, fRHS=rhs, fLHS=lhs,
              vcovOptions=from@vcovOptions, centeredVcov=from@centeredVcov,
              isEndo=from@isEndo, varNames=from@varNames,omit=from@omit,
              survOptions=from@survOptions, sSpec=from@sSpec, smooth=from@smooth)
      })

setAs("linearModel", "functionModel",
      function(from) {
          spec <- modelDims(from)          
          x <- from
          theta0 <- rep(0,spec$k)
          names(theta0) <- spec$parNames
          fct <- function(theta, x)
              {
                  names(theta) <- modelDims(x)$parNames
                  gt <- evalMoment(x, theta)
              }
          dfct <- function(theta, x)
              {
                  names(theta) <- modelDims(x)$parNames
                  gt <- evalDMoment(x, theta)
              }
          new("functionModel", X=x, fct=fct, dfct=dfct,  vcov=from@vcov,
              theta0=theta0, n=spec$n, q=spec$q, k=spec$k, parNames=names(theta0),
              momNames=spec$momNames,vcovOptions=from@vcovOptions,
              centeredVcov=from@centeredVcov,omit=integer(),survOptions=from@survOptions,
              sSpec=from@sSpec, smooth=from@smooth)
      })

setAs("allNLModel", "functionModel",
      function(from) {
          spec <- modelDims(from)          
          x <- from
          fct <- function(theta, x)
              {
                  names(theta) <- modelDims(x)$parNames
                  gt <- evalMoment(x, theta)
              }
          dfct <- function(theta, x)
              {
                  names(theta) <- modelDims(x)$parNames
                  gt <- evalDMoment(x, theta)
              }
          new("functionModel", X=x, fct=fct, dfct=dfct,  vcov=from@vcov,
              theta0=from@theta0, n=spec$n, q=spec$q, k=spec$k,
              parNames=names(from@theta0),
              momNames=spec$momNames, vcovOptions=from@vcovOptions,
              centeredVcov=from@centeredVcov,omit=integer(),
              survOptions=from@survOptions, sSpec=from@sSpec, smooth=from@smooth)
      })

setAs("slinearModel", "linearModel",
      function(from) {
          spec <- modelDims(from)
          eqnNames <- from@eqnNames
          neqn <- length(eqnNames)
          datX <- lapply(1:neqn,
                         function(i) {
                             v <- from@varNames[[i]]
                             chk <- "(Intercept)" %in% v
                             v <- v[v!="(Intercept)"]
                             X <- from@data[,v, drop=FALSE]
                             colnames(X) <- paste(eqnNames[[i]],".", v, sep="")
                             if (chk)
                                 {
                                  X <- cbind(1, X)
                                  colnames(X)[1]<-paste(eqnNames[[i]], ".Intercept", sep="")
                                 }
                             X})
          datZ <- lapply(1:neqn,
                         function(i) {
                             v <- all.vars(from@instT[[i]])
                             chk <- attr(from@instT[[i]], "intercept")==1
                             Z <- from@data[,v, drop=FALSE]
                             colnames(Z) <- paste(eqnNames[[i]],".", v, sep="")
                             if (chk)
                                 {
                                  Z <- cbind(1, Z)
                                  colnames(Z)[1]<-paste(eqnNames[[i]], ".Intercept", sep="")
                                 }
                             Z})
          nZ <- do.call("c", lapply(datZ, colnames))
          nX <- do.call("c", lapply(datX, colnames))
          datZ <- .GListToMat(datZ)
          datX <- .GListToMat(datX)
          Y <- do.call("c", modelResponse(from))
          colnames(datZ) <- nZ
          colnames(datX) <- nX
          dat <- cbind(Y, datZ, datX)
          dat <- dat[,unique(colnames(dat))]
          dat <- data.frame(dat, row.names=1:nrow(datZ))
          g <- paste("Y~", paste(nX, collapse="+"), "-1")
          g <- formula(g, .GlobalEnv)
          h <- paste("~", paste(nZ, collapse="+"), "-1")
          h <- formula(h, .GlobalEnv)
          res <- momentModel(g, h, vcov=from@vcov, vcovOptions=from@vcovOptions,
                          centeredVcov=from@centeredVcov, data=dat)
      })

setAs("slinearModel", "snonlinearModel",
      function(from) {
          spec <- modelDims(from)
          X <- model.matrix(from)
          theta0 <- rep(1,sum(spec$k))         
          names(theta0) <- paste("theta", 1:sum(spec$k), sep="")
          eqNames <- paste("Eqn", 1:length(X), sep="")
          xk <- c(0,cumsum(from@k))
          theta0 <- lapply(1:length(X), function(i) theta0[(1+xk[i]):(xk[i+1])])
          parNames <- lapply(theta0, names)
          rhs <- lapply(1:length(X), function(i){
              n <- paste("*", colnames(X[[i]]), sep="")
              n[n=="*(Intercept)"] <- ""
              n <- paste(names(theta0[[i]]), n, sep="")
              parse(text=paste(n, collapse="+", sep=""))
              })
          lhs <- lapply(1:length(X), function(i)
              parse(text=from@modelT[[i]][[2]]))
          varNames <- lapply(1:length(lhs), function(i) {
              v1 <- all.vars(lhs[[i]])
              v2 <- all.vars(rhs[[i]])
              v2 <- v2[!(v2%in%names(theta0[[i]]))]
              c(v1,v2)})
              
          Y <- do.call(cbind, modelResponse(from))
          colnames(Y) <- sapply(lhs, all.vars)
          X <- do.call(cbind, X)
          X <- X[,!duplicated(colnames(X))]
          X <- X[,colnames(X)!="(Intercept)"]
          Z <- do.call(cbind, model.matrix(from, type="instruments"))
          Z <- Z[,!duplicated(colnames(Z))]
          Z <- Z[,colnames(Z) != "(Intercept)"]
          dat <- cbind(X, Y[,!(colnames(Y) %in% colnames(X))])
          dat <- cbind(dat, Z[,!(colnames(Z)%in%colnames(dat))])
          new("snonlinearModel", data=as.data.frame(dat), instT=from@instT,
              vcov=from@vcov, theta0=theta0, n=spec$n, q=spec$q,k=spec$k,
              parNames=parNames, momNames=from@momNames, fRHS=rhs,
              fLHS=lhs, eqnNames=eqNames, vcovOptions=from@vcovOptions,
              centeredVcov=from@centeredVcov, sameMom=from@sameMom,
              SUR=from@SUR, varNames=varNames, isEndo=from@isEndo,
              omit=from@omit, survOptions=from@survOptions,
              sSpec=from@sSpec, smooth=from@smooth)
      })


setAs("sysMomentWeights", "momentWeights",
      function(from) {
          w <- quadra(from)
          if (is.character(w))
              w <- "ident"
          new("momentWeights", w=w, type="weights", wSpec=list())
      })
          

setAs("rslinearModel", "rlinearModel", 
     function(from) {
          m <- as(from, "slinearModel")
          m <- as(m, "linearModel")
          restModel(m, from@cstLHS, from@cstRHS)
      })
