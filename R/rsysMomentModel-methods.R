######## Methods for restricted System of Equations
######################################################

## Hidden function to arrange the restrictions

.imposeSNLRestrict <- function(R, object)
{
    spec <- modelDims(object)
    allPar <- do.call("c", c(list(),spec$parNames))
    theta0 <- do.call("c", c(list(),spec$theta0))
    fRHS <- object@fRHS
    fLHS <- object@fLHS
    chk <- sapply(R, function(r) all(all.vars(r) %in% allPar))
    chk2 <- sapply(R, function(r) sum(sapply(spec$parNames,
                                             function(ni) any(all.vars(r) %in% ni))))
    crossRest <- any(chk2>1)
    if (!all(chk))
        stop("Wrong coefficient names in some of the restrictions")
    rest <- sapply(R, function(r) as.character(r[[2]]))
    if (!all(sapply(rest, function(x) length(x)==1)))
        stop("LHS of R formulas must contain only one coefficient")
    dR <-numeric()
    for (r in R)
    {
        lhs <- sapply(allPar, function(pn)
            eval(D(r[[2]], pn), as.list(theta0)))
        rhs <- sapply(allPar, function(pn)
            eval(D(r[[3]], pn), as.list(theta0)))
        dR <- rbind(dR, lhs-rhs)
    }
    if (any(is.na(dR)) || any(!is.finite(dR)))
        stop("The derivative of the constraints at theta0 is either infinite or NAN")
    if (qr(dR)$rank < length(R))
        stop("The matrix of derivatives of the constraints is not full rank")
    rhs <- object@fRHS
    lhs <- object@fLHS
    env <-  new.env()
    varNames <- list()
    for (r in R)
    {
        assign(as.character(r[2]),
               parse(text=paste("(", as.character(r[3]), ")", sep=""))[[1]], env)
    }
    for (i in 1:length(rhs))
    {
        rhs[[i]] <- as.expression(do.call('substitute', list(rhs[[i]][[1]], env)))
        if (!is.null(lhs[[i]]))
            lhs[[i]] <- as.expression(do.call('substitute', list(lhs[[i]][[1]], env)))
        tmp <- c(all.vars(rhs[[i]]))
    }
    k <- sum(spec$k)-length(R)
    parNames <- if (is.list(spec$parNames))
                    {
                        lapply(spec$parNames, function(pi) pi[!(pi %in% rest)])
                    } else {
                        spec$parNames[!(spec$parNames %in% rest)]
                    }
    theta0 <- if (is.list(spec$theta0))
              {
                  lapply(spec$theta0, function(ti) ti[!(names(ti) %in% rest)])
              } else {
                  spec$theta0[!(names(spec$theta0) %in% rest)]
              }
    
    if (length(rhs) == 1)
    {
        rhs <- rhs[[1]]
        lhs <- lhs[[1]]
    }
    if (is.list(parNames))
        k <- sapply(parNames, length)
    list(rhs=rhs, lhs=lhs, parNames=parNames,
         theta0=theta0, k=k, crossEquRest=crossRest)
}


## Constructor of restricted models

setMethod("restModel", signature("slinearModel"),
          function(object, R, rhs=NULL)
          {
              parNames <- paste(rep(object@eqnNames, object@k), ".",
                                do.call("c", object@parNames), sep = "")
              
              if (is.character(R))
              {
                  parNames <- gsub("\\(Intercept\\)", "Intercept", parNames)
                  R <- gsub("\\(Intercept\\)", "Intercept", R)
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


setMethod("restModel", signature("snonlinearModel"),
          function(object, R, rhs=NULL)
          {
              if (!is.null(rhs))
                  warning("rhs is ignored for nonlinear models")
              if (is.character(R))
                  {
                      R2 <- list()
                      R <- gsub("=", "~", R, fixed=TRUE)
                      for (r in R)
                          R2 <- c(R2, as.formula(r, .GlobalEnv))
                      R <- R2
                  } else {
                      if (!is.list(R))
                          {
                              if(!inherits(R,"formula"))
                                  stop("R must be a formula or a list of formulas")
                              R <- list(R)
                          } else {
                              chk <- sapply(R, function(r) inherits(r,"formula"))
                              if (!all(chk))
                                  stop("R must be a formula, a list of formulas or a vector of characters")
                          }
                  }
              res <- .imposeSNLRestrict(R, object)
              cstSpec <- list(newParNames = res$parNames,
                              originParNames=object@parNames,
                              k=res$k, theta0=res$theta0, fRHS=res$rhs, fLHS=res$lhs,
                              crossEquRest=res$crossEquRest)
              new("rsnonlinearModel",  R=R, cstSpec=cstSpec, object)
          })

## coef

setMethod("coef","rslinearModel",
          function (object, theta) 
          {
              spec <- modelDims(object)
              theta <- setCoef(object, theta)
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

setMethod("coef", "rsnonlinearModel",
          function(object, theta)
          {
              spec <- modelDims(object)
              theta <- setCoef(object, theta)
              names(theta) <- NULL
              theta <- do.call("c", theta)
              theta2 <- rep(0,sum(object@k))
              names(theta2) <- do.call("c", object@parNames)
              theta2[names(theta)] <- theta
              chk <- sapply(object@R, function(r) is.numeric(r[[3]]))
              for (r in object@R[chk])
                  theta2[as.character(r[[2]])] <- r[[3]]
              for (r in object@R[!chk])
                  theta2[as.character(r[[2]])] <- eval(r[[3]], as.list(theta2))
              .tetReshape(theta2, spec$eqnNames, object@parNames)
          })

## modelDims

setMethod("modelDims", "rslinearModel",
          function(object) {
              res <- object@cstSpec
              list(k = res$k, q = res$q, n = res$n, parNames = res$newParNames, 
                   momNames = object@momNames, eqnNames=res$eqnNames,
                   isEndo=res$isEndo)
          })

setMethod("modelDims", "rsnonlinearModel",
          function(object) {
              cst <- object@cstSpec
              k <- cst$k
              parNames <- cst$newParNames
              theta0 <- cst$theta0
              fRHS <- cst$fRHS
              fLHS <- cst$fLHS
              list(k=k, q=object@q, n=object@n, parNames=parNames,
                   momNames=object@momNames, theta0=theta0,
                   fRHS=fRHS, fLHS=fLHS, eqnNames=object@eqnNames,
                   isEndo=object@isEndo)
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

setMethod("printRestrict", "rsnonlinearModel",
          function(object){
              cat("Constraints:\n")
              parNames <- object@parNames
              eqn <- object@eqnNames
              chk <- lapply(object@R, function(ri)
                  which(sapply(parNames, function(pi) any(all.vars(ri) %in% pi))))
              for (i in unique(chk))
              {
                  chk2 <- which(sapply(chk, function(ci) isTRUE(all.equal(ci,i)))) 
                  nr <- length(i)
                  if (nr>1)
                  {
                      cat("\tInvolving equations: ")
                      if (nr>2)
                          sep <- c(rep(", ", nr-2), " and ", "")
                      else
                          sep <- c(" and ", "")                      
                      cat(paste(eqn[i], sep, collapse="", sep=""), "\n")
                  } else {
                      cat("\tInvolving equation: ", eqn[i], "\n", sep="")
                  }
                  for (ri in object@R[chk2])
                      {
                          cat("\t\t")
                          print(ri)
                      }
              }
              cat("\n")
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

setMethod("print", "rsnonlinearModel", 
          function (x, ...) {
              callNextMethod()
              cat("(The number of endogenous variables is unreliable)\n")              
              printRestrict(x)
              cat("\n")
          } )


## getRestrict

setMethod("getRestrict", "sysModel",
          function(object, theta, R, rhs=NULL) {
              robject <- restModel(object, R, rhs)
              getRestrict(robject, theta)
          })

setMethod("getRestrict", "rslinearModel",
          function(object, theta) {
              theta <- setCoef(as(object, "slinearModel"), theta)
              theta <- do.call("c", theta)
              R <- c(object@cstLHS%*%theta)
              parNames <- paste(rep(object@eqnNames, object@k), ".",
                                do.call("c", object@parNames), sep = "")
              cst <- .printHypothesis(object@cstLHS, object@cstRHS, parNames)
              list(dR=object@cstLHS, R=R, q=object@cstRHS, hypo=cst,
                   orig.R=object@cstLHS, orig.rhs=object@cstRHS)
          })

setMethod("getRestrict", "rsnonlinearModel",
          function(object, theta) {
              theta <- setCoef(as(object, "snonlinearModel"), theta)              
              names(theta) <- NULL
              theta <- do.call("c", theta)
              dR <-numeric()
              R <- numeric()
              parNames <- names(theta)
              for (r in object@R)
                  {
                      dlhs <- sapply(parNames, function(pn)
                          eval(D(r[[2]], pn), as.list(theta)))
                      drhs <- sapply(parNames, function(pn)
                          eval(D(r[[3]], pn), as.list(theta)))
                      dR <- rbind(dR, dlhs-drhs)
                      lhs <- eval(r[[2]], as.list(theta))
                      rhs <- eval(r[[3]], as.list(theta))
                      R <- c(R, lhs-rhs)                      
                  }
              if (any(is.na(c(R,dR))) || any(!is.finite(c(dR,R))))
                  stop("Some values in R or dR at theta are either infinite or NAN")
              if (qr(dR)$rank < length(R))
                  stop("The matrix of derivatives of the constraints is not full rank")
              hypo <- sapply(object@R, function(r) capture.output(print(r)))
              list(dR=dR, R=R, q=rep(0, nrow(dR)), hypo=hypo,
                   orig.R=object@R, orig.rhs=NULL)
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


setMethod("[", c("rsnonlinearModel", "numeric", "missing"),
          function(x, i, j)
          {
              i <- as.integer(i)
              if (x@cstSpec$crossEquRest)
                  stop("Cannot select an equation when the system has cross-equation restrictions")
              if (length(i) > 1)
              {
                  warning("You can only select one equation, the first element of i is used")
                  i <- i[1]
              }
              m <- as(x, "snonlinearModel")[i]
              chk <- sapply(x@R, function(ri) any(all.vars(ri) %in% m@parNames))
              R <- x@R[chk]
              if (length(R)==0)
                  return(m)
              restModel(m, R)
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

setMethod("evalMoment", "rsysModel",
          function (object, theta) 
          {
              theta <- coef(object, theta)
              evalMoment(as(object, "sysModel"), theta)
          })


## evalDMoment

setMethod("evalDMoment","rslinearModel",
          function (object, theta) 
          {
              chk <- ifelse(is.null(object@cstSpec$crossEquRest), FALSE,
                            object@cstSpec$crossEquRest)
              neqn <- length(object@eqnNames)
              if (!chk)
                  dgt <- lapply(1:neqn, function(i) evalDMoment(object[i]))
              else 
                  dgt <- list(neqn*evalDMoment(as(object, "rlinearModel")))
              names(dgt) <- modelDims(object)$eqnNames
              dgt
          })

setMethod("evalDMoment", signature("rsnonlinearModel"),
          function(object, theta, impProb=NULL)
          {
              De <- Dresiduals(object, theta)
              Z <- model.matrix(as(object, "snonlinearModel"), "instrument")
              spec <- modelDims(object)
              if (is.null(impProb))
                  impProb <- 1/spec$n
              sG <- lapply(1:length(spec$eqnNames), function(eq) {
                  G <- apply(De[[eq]],2, function(x)
                  {
                      tmp <- Z[[eq]]*x
                      if (object@smooth)
                          tmp <- stats::kernapply(tmp, object@sSpec@w)
                      colSums(tmp*impProb)
                  })
                  G <- as.matrix(G)
                  dimnames(G) <- list(spec$momNames[[eq]], colnames(De[[eq]]))
                  G
              })
              names(sG) <- spec$eqnNames
              sG
          })

## residuals

setMethod("residuals", "rsysModel",
          function (object, theta) 
          {
              theta <- coef(object, theta)
              residuals(as(object, "sysModel"), theta)
          })

## model.matrix

setMethod("model.matrix", "rsnonlinearModel",
          function(object, type =  c("regressors", "instruments"))
          {
              type <- match.arg(type)
              model.matrix(as(object, "snonlinearModel"), type)
          })


## Dresiduals

setMethod("Dresiduals", signature("rsnonlinearModel"),
          function(object, theta) {
              spec <- modelDims(object)
              theta <- setCoef(object, theta)
              if (!object@cstSpec$crossEquRest)
              {
                  De <- Dresiduals(as(object, "snonlinearModel"), coef(object, theta))
                  for (i in 1:length(De))
                      De[[i]] <- De[[i]][,colnames(De[[i]]) %in% spec$parNames[[i]]]
                  return(De)
              }
              neqn <- length(object@eqnNames)
              fLHS <- spec$fLHS
              fRHS <- spec$fRHS
              names(theta) <- NULL
              theta <- do.call("c", theta)
              nt <- names(theta)
              varList <- c(as.list(theta), as.list(object@data))
              sDe <- lapply(1:neqn, function(eq){
                  De <- numeric()
                  for (i in nt)
                  {
                      if (!is.null(fLHS[[eq]]))
                          d <- eval(D(fLHS[[eq]], i), varList)      
                      else
                          d <- 0
                      De <-  cbind(De, d-matrix(eval(D(fRHS[[eq]], i), varList),spec$n,1))
                  }
                  colnames(De) <- nt
                  De
              })
              names(sDe) <- spec$eqnNames
              sDe
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

