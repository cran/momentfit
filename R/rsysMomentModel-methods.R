######## Methods for restricted System of Equations
######################################################

## Constructor of restricted models

setMethod("restModel", signature("slinearModel"),
          function(object, R, rhs=NULL)
          {
              parNames <- paste(rep(object@eqnNames, object@k), ".",
                                do.call("c", object@parNames), sep = "")                 
              if (is.character(R))
              {
                  res <- .makeHypothesis(parNames, R, rhs)
                  R <- res$R
                  rhs <- res$rhs
              } else if (is.list(R)) {
                  nrest <- sum(sapply(R, length))
                  R2 <- numeric()
                  rhs <- numeric()
                  done <- 0
                  for (i in 1:length(R))
                  {
                      if (length(R[[i]]) == 0)
                      {
                          R2 <- cbind(R2, matrix(0, nrest, object@k[i]))
                      } else {
                          tmp <- .makeHypothesis(object@parNames[[i]], R[[i]], NULL)
                          rhs <- c(rhs, tmp$rhs)
                          tmp2 <- matrix(0, nrest, object@k[i])
                          tmp2[(done+1):(done+length(R[[i]])),] <- tmp$R
                          done <- done+length(R[[i]])
                          R2 <- cbind(R2, tmp2)
                      }
                  }
                  R <- R2
              } else {
                  if (is.null(rhs))
                      rhs <- rep(0,nrow(R))
              }
              neqn <- length(object@eqnNames)
              ind <- rbind(c(1, cumsum(object@k)[-neqn]+1),
                           cumsum(object@k))
              chk <- sapply(1:nrow(R), function(i)
                  sapply(1:neqn, function(j) any(R[i,seq(ind[1,j],ind[2,j])]!=0)))
              res <- try(.imposeRestrict(R,rhs,parNames), silent=TRUE)
              if (any(class(res) == "try-error"))
                  stop("Failed to construct restricted model from the provided restrictions can you simplify it?")
              res$crossEquRest <- any(colSums(chk)>1)
              res$eqnSelect <- ind
              isEndo <- object@isEndo
              rtet <- res$theta
              isEndo <- c(crossprod(do.call("c", isEndo), rtet)) != 0
              if (!res$crossEquRest)
              {
                  res$k <- sapply(object@eqnNames,
                                  function(en) length(grep(en, res$newParNames)))
                  res$newParNames <- lapply(object@eqnNames, function(en)
                      gsub(paste(en,".",sep=""), "",
                           res$newParNames[grep(en, res$newParNames)]))
                  res$n <- object@n
                  res$q <- object@q
                  res$eqnNames <- object@eqnNames
                  res$isEndo <- .tetReshape(isEndo, res$eqnNames, res$newParNames)
              } else {
                  res$n <- object@n*length(object@eqnNames)
                  res$q <- sum(object@q)
                  res$eqnNames <- "combinedEqns"
                  names(isEndo) <- res$newParNames
                  res$isEndo <- list(isEndo)
                  res$newParNames <- list(res$newParNames)
                  names(res$isEndo) <- names(res$newParNames) <- res$eqnNames 
              }
              res$originParNames <- object@parNames
              new("rslinearModel",  cstLHS=R, cstRHS=rhs,
                  cstSpec=res, object)
          })

## coef

setMethod("coef","rslinearModel",
          function (object, theta) 
          {                  
              spec <- modelDims(object)
              if (!is.list(theta))
                  stop("theta must be a list")
              if (length(theta) != length(spec$eqnNames)) 
                  stop("Wrong number of coefficients")
              if (!object@cstSpec$crossEquRest)
              {
                  tet <- lapply(1:length(theta), function(i)
                      coef(object[i], theta[[i]]))
                  names(tet) <- object@eqnNames
                  return(tet)
              }
              theta <- coef(as(object, "rlinearModel"), theta[[1]])
              theta <- .tetReshape(theta, object@eqnNames, object@parNames)
              theta
          })

## modelDims

setMethod("modelDims", "rslinearModel",
          function(object) {
              res <- object@cstSpec
              list(k = res$k, q = res$q, n = res$n, parNames = res$newParNames, 
                   momNames = object@momNames, eqnNames=res$eqnNames,
                   isEndo=res$isEndo)
          })

## printRestrict

setMethod("printRestrict", "rslinearModel",
          function(object){
              parNames <- paste(rep(object@eqnNames, object@k), ".",
                                do.call("c", object@parNames), sep = "")
              cst <- .printHypothesis(object@cstLHS, object@cstRHS, parNames)
              cat("Constraints:\n")
              for (i in 1:length(cst))
                  cat("\t", cst[i], "\n")
          })

## print

setMethod("print", "rslinearModel", 
          function (x, ...) {
              callNextMethod()
              if (!x@cstSpec$crossEquRest)
              {
                  cat("**Equation by Equation restrictions**\n")
                  for (i in 1:length(x@eqnNames))
                  {
                      m <- x[i]
                      if (inherits(m, "rlinearModel"))
                      {
                          cat("**", x@eqnNames[i], "**\n", sep="")
                          printRestrict(m)
                          cat("\n")
                      }
                  }
              } else {
                  printRestrict(x)
              }
          })

## getRestrict

setMethod("getRestrict", "sysModel",
          function(object, theta, R, rhs=NULL) {
              robject <- restModel(object, R, rhs)
              getRestrict(robject, theta)
          })

setMethod("getRestrict", "rslinearModel",
          function(object, theta) {
              theta <- do.call("c", theta)
              R <- c(object@cstLHS%*%theta)
              parNames <- paste(rep(object@eqnNames, object@k), ".",
                                do.call("c", object@parNames), sep = "")
              cst <- .printHypothesis(object@cstLHS, object@cstRHS, parNames)
              list(dR=object@cstLHS, R=R, q=object@cstRHS, hypo=cst,
                   orig.R=object@cstLHS, orig.rhs=object@cstRHS)
          })

### Subset

setMethod("[", c("rslinearModel", "numeric", "missing"),
          function(x, i, j)
          {
              i <- as.integer(i)
              if (length(i) > 1)
              {
                  warning("You can only select one equation, the first element of i is used")
                  i <- i[1]
              }
              if (x@cstSpec$crossEquRest)
              {
                  if (i !=1L)
                      stop("There is only one equation with cross-equation restrictions")
                  return(as(x, "rlinearModel"))
              }
              m <- as(x, "slinearModel")[i]
              sel <- x@cstSpec$eqnSelect
              R <- x@cstLHS[, sel[1,i]:sel[2,i]]
              chk <- apply(R, 1, function(x) all(x==0))
              if (all(chk))
                  return(m)
              R <- R[!chk,,drop=FALSE]
              rhs <- x@cstRHS[!chk]
              restModel(m, R, rhs)
          })

### model.matrix

setMethod("model.matrix", "rslinearModel",
          function (object, type = c("regressors", "instruments")) 
          {
              type <- match.arg(type)
              if (!object@cstSpec$crossEquRest)
                  return(callNextMethod())
              if (type == "instruments")
                  mm <- model.matrix(as(object, "slinearModel"), type="instruments")
              else
                  mm <- list(model.matrix(as(object, "rlinearModel")))
              names(mm) <- modelDims(object)$eqnNames
              mm
          })

## modelResponse

setMethod("modelResponse", "rslinearModel",
          function(object) {
              if (!object@cstSpec$crossEquRest)
                  return(callNextMethod())
              mr <- list(modelResponse(as(object, "rlinearModel")))
              names(mr) <- modelDims(object)$eqnNames
              mr
          })

## evalMoment

setMethod("evalMoment", "rslinearModel",
          function (object, theta) 
          {
              theta <- coef(object, theta)
              evalMoment(as(object, "slinearModel"), theta)
          })

## evalDMoment

setMethod("evalDMoment","rslinearModel",
          function (object, theta) 
          {
              neqn <- length(object@eqnNames)
              if (!object@cstSpec$crossEquRest)
                  dgt <- lapply(1:neqn, function(i) evalDMoment(object[i]))
              else 
                  dgt <- list(neqn*evalDMoment(as(object, "rlinearModel")))
              names(dgt) <- modelDims(object)$eqnNames
              dgt
          })

## residuals

setMethod("residuals", "rslinearModel",
          function (object, theta) 
          {
              theta <- coef(object, theta)
              residuals(as(object, "slinearModel"), theta)
          })

## solveGmm

setMethod("solveGmm", c("rslinearModel","sysMomentWeights"),
          function (object, wObj, theta0 = NULL) 
          {
              if (object@cstSpec$crossEquRest)
              {
                  res <- solveGmm(as(object, "rlinearModel"), as(wObj, "momentWeights"))
                  res$theta <- list(res$theta)
                  names(res$theta) <- modelDims(object)$eqnNames
                  return(res)
              }
              if (wObj@type == "iid" && object@sameMom) 
                  return(ThreeSLS(object, Sigma = wObj@Sigma, qrZ = wObj@w, 
                                  coefOnly = TRUE))
              spec <- modelDims(object)
              Y <- modelResponse(object)
              Z <- model.matrix(object, type = "instruments")
              G <- evalDMoment(object)
              Syz <- lapply(1:length(Y), function(i) colMeans(Y[[i]]*Z[[i]]))
              Syz <- do.call("c", Syz)
              G <- momentfit:::.GListToMat(G)
              T1 <- quadra(wObj, G)
              T2 <- quadra(wObj, G, Syz)
              theta <- -solve(T1, T2)
              spec <- modelDims(object)
              theta <- momentfit:::.tetReshape(theta, object@eqnNames, spec$parNames)
              list(theta = theta, convergence = NULL)
          })


### Three Stage Least Squares

setMethod("ThreeSLS","rslinearModel",
          function(model, coefOnly=FALSE, qrZ=NULL, Sigma=NULL) {
              if (!model@cstSpec$crossEquRest)
                  callNextMethod()
              else
                  stop("Systems with cross-equation restrictons are not considered system of equations and cannot be estimated by 3SLS")
          })

## evalWeights, we just need to transform the model with cross-equation restrictions
## back to a system of equation first.

setMethod("evalWeights", "rslinearModel",
          function(object, theta = NULL, w="optimal", wObj=NULL)
          {
              if (object@cstSpec$crossEquRest)
              {
                  if (!is.null(theta))
                      theta <- coef(object, theta)
                  return(evalWeights(as(object, "slinearModel"), theta, w, wObj))
              }
              callNextMethod()
          })

## gmmFit. Almost like sysModes method, bu we need to check a few things

setMethod("gmmFit", signature("rslinearModel"), valueClass="sgmmfit", 
          function(model, type=c("twostep", "iter","cue", "onestep"),
                   itertol=1e-7, initW=c("ident", "tsls", "EbyE"),
                   weights="optimal", itermaxit=100,
                   efficientWeights=FALSE, theta0=NULL, EbyE=FALSE, ...)
          {
              type <- match.arg(type)
              initW <- match.arg(initW)
              if (model@cstSpec$crossEquRest)
              {
                  if (EbyE || initW=="EbyE")
                      stop("Models with cross-equation restrictions cannot be estimated equation by equation because it is not considered a system of equations")
                  if (initW == "tsls")
                      stop("Models with cross-equation restrictions cannot be initiated by an equation by equation 2SLS because it is not considered a system of equations")
                  if (type=="onestep" || (is.character(weights) && weights=="ident"))
                  {
                      wObj <- evalWeights(model, w="ident")
                      return(gmmFit(model, w=wObj))
                  }
              }
              callNextMethod()
          })

