### A function to help imposing the constraint
### the new Y is Y-c(X%*%minY)
### and the new X is X%*%theta

.getRestNames <- function(theta, parNames)
{
    X <- character()
    if (ncol(theta) == 0)
        return(X)
    for (i in 1:ncol(theta))
    {
        t <- theta[w <- which(theta[,i]!=0),i]
        tmp <- c(as.character(t[1]))
        if (length(t)>1)
            tmp <- c(tmp, ifelse(t[-1]<0, paste("-", abs(t[-1]), sep=""),
                                 paste("+", t[-1], sep="")))
        if (any(abs(t) == 1))
            tmp[abs(t)==1] <- gsub("1","",tmp[abs(t)==1])
        X[i] <- paste(tmp, parNames[w], collapse="", sep="")
        if (length(t)>1)
            {
                X[i] <- paste("(",X[i],")",sep="")
            }else if (t!=1) {
                X[i] <- paste("(",X[i],")",sep="")
            }
    }
    X
}

.imposeRestrict <- function(R,q,parNames)
{
    chk  <- apply(R, 1, function(x) sum(x!=0)==1)
    minY <- rep(0,ncol(R))
    theta <- diag(ncol(R))
    done <- which(chk)
    while(TRUE)
    {
        if (all(!chk))
            break
        r <- R[chk,,drop=FALSE]
        ri <- apply(r,1,function(x) which(x!=0))
        b <- q[chk]/r[cbind(1:nrow(r),ri)]
        minY[ri] <- b
        diag(theta)[ri] <- 0
        w <- which(!chk)
        r2 <- R[w,ri,drop=FALSE]
        q[w] <- q[w] - c(r2%*%b)
        R[w,ri] <- 0
        R[chk,] <- 0
        chk  <- apply(R, 1, function(x) sum(x!=0)==1)
        done <- c(done, which(chk))
    }
    if (length(done) == nrow(R))
    {
        theta <- theta[,apply(theta,2,function(x) any(x!=0)), drop=FALSE]
        newParNames <- .getRestNames(theta, parNames)
        return(list(theta=theta, minY=minY, newParNames=newParNames,
             originParNames=parNames, k=length(newParNames)))
    }
    if (length(done)>0)
    {
        todo <- (1:nrow(R))[-sort(done)]
    } else {
        todo <- 1:nrow(R)
    }
    for (i in todo)
    {
        r <- R[i,]
        t1 <- which(r!=0)
        st1 <- 1
        while (sum(theta[,t1[st1]]!=0)>1)
            st1 <- st1+1
        q[i] <- q[i]/r[t1[st1]]
        r <- r/r[t1[st1]]
        diag(theta)[t1[st1]] <- 0
        minY[t1[st1]] <- q[i]
        theta[t1[st1],t1[-st1]] <- -r[t1[-st1]]
    }
    theta <- theta[,apply(theta,2,function(x) any(x!=0)), drop=FALSE]
    newParNames <- .getRestNames(theta, parNames)
    list(theta=theta, minY=minY, newParNames=newParNames,
         originParNames=parNames, k=length(newParNames))
}

.imposeNLRestrict <- function(R, object)
    {
        chk <- sapply(R, function(r) all(all.vars(r) %in% object@parNames))
        if (!all(chk))
            stop("Wrong coefficient names in some of the restrictions")
        rest <- sapply(R, function(r) as.character(r[[2]]))
        if (!all(sapply(rest, function(x) length(x)==1)))
            stop("LHS of R formulas must contain only one coefficient")
        dR <-numeric()
        for (r in R)
            {
                lhs <- sapply(object@parNames, function(pn)
                    eval(D(r[[2]], pn), as.list(object@theta0)))
                rhs <- sapply(object@parNames, function(pn)
                    eval(D(r[[3]], pn), as.list(object@theta0)))
                dR <- rbind(dR, lhs-rhs)
            }
        if (any(is.na(dR)) || any(!is.finite(dR)))
            stop("The derivative of the constraints at theta0 is either infinite or NAN")
        if (qr(dR)$rank < length(R))
            stop("The matrix of derivatives of the constraints is not full rank")
        rhs <- object@fRHS
        lhs <- object@fLHS
        env <-  new.env()        
        for (r in R)
        {
            assign(as.character(r[2]),
                   parse(text=paste("(", as.character(r[3]), ")", sep=""))[[1]], env)
        }
        rhs <- as.expression(do.call('substitute', list(rhs[[1]], env)))
        if (!is.null(lhs))
            lhs <- as.expression(do.call('substitute', list(lhs[[1]], env)))
        k <- object@k-length(R)
        parNames <- object@parNames[!(object@parNames %in% rest)]
        theta0 <- object@theta0[!(object@parNames %in% rest)]        
        list(rhs=rhs, lhs=lhs, parNames=parNames, theta0=theta0, k=k)
    }

.imposefRestrict <- function(R, object)
    {
        chk <- sapply(R, function(r) all(all.vars(r) %in% object@parNames))
        if (!all(chk))
            stop("Wrong coefficient names in some of the restrictions")
        rest <- sapply(R, function(r) as.character(r[[2]]))
        if (!all(sapply(rest, function(x) length(x)==1)))
            stop("LHS of R formulas must contain only one coefficient")
        k <- object@k-length(R)
        parNames <- object@parNames[!(object@parNames %in% rest)]
        theta0 <- object@theta0[!(object@parNames %in% rest)]        
        list(parNames=parNames, theta0=theta0, k=k)
    }

.imposeFORMRestrict <- function(R, object)
    {
        chk <- sapply(R, function(r) all(all.vars(r) %in% object@parNames))
        if (!all(chk))
            stop("Wrong coefficient names in some of the restrictions")
        rest <- sapply(R, function(r) as.character(r[[2]]))
        if (any(duplicated(rest)))
            stop("LHS of R must not have duplicated variables")
        if (!all(sapply(rest, function(x) length(x)==1)))
            stop("LHS of R formulas must contain only one coefficient")
        dR <-numeric()
        for (r in R)
            {
                lhs <- sapply(object@parNames, function(pn)
                    eval(D(r[[2]], pn), as.list(object@theta0)))
                rhs <- sapply(object@parNames, function(pn)
                    eval(D(r[[3]], pn), as.list(object@theta0)))
                dR <- rbind(dR, lhs-rhs)
            }
        if (any(is.na(dR)) || any(!is.finite(dR)))
            stop("The derivative of the constraints at theta0 is either infinite or NAN")
        if (qr(dR)$rank < length(R))
            stop("The matrix of derivatives of the constraints is not full rank")
        rhs <- object@fRHS
        lhs <- object@fLHS
        env <-  new.env()
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
        }
        k <- object@k-length(R)
        parNames <- object@parNames[!(object@parNames %in% rest)]
        if (length(parNames)!=k)
            stop("Failed to create the restricted model")
        theta0 <- object@theta0[!(object@parNames %in% rest)]        
        list(rhs=rhs, lhs=lhs, parNames=parNames, theta0=theta0, k=k)
    }

################## model.matrix and modelResponse #################
### I did not make model.response as generic because it is not
### a method in stats and I want different arguments

setMethod("modelResponse", signature("rlinearModel"),
          function(object)
          {
              Y <- model.response(object@modelF)
              minY <- object@cstSpec$minY
              if (all(minY==0))
                  return(Y)
              ti <- attr(object@modelF, "terms")
              X <- model.matrix(ti, object@modelF)[,minY!=0]
              minY <- minY[minY!=0]
              Y <- Y - colSums(t(X)*minY)
              Y
          })

setMethod("model.matrix", signature("rlinearModel"),
          function(object, type=c("regressors","instruments"))
          {
              type <- match.arg(type)
              if (type == "instruments")
              {
                  mat <- callNextMethod(object, type=type)
              } else {
                  res <- object@cstSpec
                  theta <- res$theta
                  ti <- attr(object@modelF, "terms")
                  mat <- model.matrix(ti, object@modelF)[,]
                  w <- apply(theta,2,function(x) all(x==0))
                  theta <- theta[,!w,drop=FALSE]
                  mat <- mat%*%theta
                  colnames(mat) <- res$newParNames
              }
              mat
          })

############### modelDims #######################

setMethod("modelDims", "rlinearModel",
          function(object) {
              res <- object@cstSpec
              list(k=res$k, q=object@q, n=object@n, parNames=res$newParNames,
                   momNames=object@momNames, isEndo=res$isEndo)
          })

setMethod("modelDims", "rformulaModel",
          function(object) {
              res <- object@cstSpec
              list(k=res$k, q=object@q, n=object@n, parNames=res$newParNames,
                   momNames=object@momNames, theta0=res$theta0,
                   fRHS=res$fRHS, fLHS=res$fLHS)
          })

setMethod("modelDims", "rnonlinearModel",
          function(object) {
              res <- object@cstSpec
              list(k=res$k, q=object@q, n=object@n, parNames=res$newParNames,
                   momNames=object@momNames, theta0=res$theta0,
                   fRHS=res$fRHS, fLHS=res$fLHS)
          })

setMethod("modelDims", "rfunctionModel",
          function(object) {
              res <- object@cstSpec
              list(k = res$k, q = object@q, n = object@n, parNames = res$newParNames, 
                   momNames = object@momNames, theta0 = res$theta0, 
                   fct = object@fct, dfct = object@dfct)
          })

### evalDMoment

setMethod("evalDMoment", signature("rfunctionModel"),
          function(object, theta, impProb=NULL, lambda=NULL)
          {
              G <- evalDMoment(as(object, "functionModel"),
                               coef(object, theta), impProb, lambda)
              ntheta <- modelDims(object)$parNames
              G <- G[,colnames(G) %in% ntheta, drop=FALSE]
              G
          })

setMethod("evalDMoment", signature("rformulaModel"),
          function(object, theta, impProb=NULL, lambda=NULL)
          {
              G <- evalDMoment(as(object, "formulaModel"),
                               coef(object, theta), impProb, lambda)
              ntheta <- modelDims(object)$parNames
              G <- G[,colnames(G) %in% ntheta, drop=FALSE]
              G
          })

setMethod("evalDMoment", signature("rnonlinearModel"),
          function(object, theta, impProb=NULL, lambda=NULL)
          {
              G <- evalDMoment(as(object, "nonlinearModel"),
                               coef(object, theta), impProb, lambda)
              ntheta <- modelDims(object)$parNames
              G <- G[,colnames(G) %in% ntheta, drop=FALSE]
              G
          })


### print restricted equation

.printRFct <- function(object)
{
    res <- object@cstSpec
    parNames <- object@parNames
    y <- colnames(object@modelF)[1L]
    w <- which(res$minY!=0)
    minY <- res$minY[w]
    minY <- ifelse(minY<0, paste("+", abs(minY),sep=""), paste("-", minY,sep=""))
    minY <- ifelse(minY=="+1", "+", minY)
    minY <- ifelse(minY=="-1", "-", minY)
    n <- paste(minY, paste(parNames[w]), collapse="",sep="")
    if (any(minY!=0))
        lhs <- paste("(",y, n, ")", sep="")
    else
        lhs <- y
    theta <- res$theta
    X <- res$newParNames
    rhs <- paste(X, collapse="+",sep="")
    paste(lhs,"=",rhs)
}


### Tools for setting restrictions

.printHypothesis <- function (L, rhs, cnames) 
{
    hyp <- character()
    for (i in 1:nrow(L)) {
        sel <- L[i, ] != 0
        nms <- cnames[sel]
        h <- L[i,sel]
        if (abs(h[1]) == 1)
            {
                h1 <- ifelse(h[1]<0, paste("-", nms[1], sep=""), nms[1])
            } else {
                h1 <- ifelse(h[1]<0, paste("-", -h[1], sep=""),
                             as.character(h[1]))
                h1 <- paste(h1, nms[1])
            }
        if (length(h)>1)
        {
            h2 <- ifelse(h[-1] < 0,
                         paste(" -", -h[-1]),
                         paste(" +",  h[-1]))
            if (any(abs(h[-1])==1))
                h2[abs(h[-1])==1] <- gsub("1", "", h2[abs(h[-1])==1])
            h2 <- paste(h2, nms[-1], sep="")
            h1 <- paste(c(h1,h2), collapse="")
        }
        h1 <- paste(h1, "=", rhs[i])
        hyp[i] <- h1
    }
    hyp
}

.makeHypothesis <- function (cnames, hypothesis, rhs = NULL) 
{
    l <- list()
    n <- length(hypothesis)
    k <- length(cnames)
    # an attempt to rename all special variable names (from transformed I() e.g. or
    # interaction :. 
    newN <- paste("theta", 1:k, sep="")
    tmp <- cnames
    hasI <- grepl("I(", cnames, fixed=TRUE)
    for (w in which(hasI))
    {
        hypothesis <- gsub(tmp[w], newN[w], hypothesis, fixed=TRUE)
        cnames <- gsub(tmp[w], newN[w], cnames, fixed=TRUE)
    }
    for (w in which(!hasI))
        {
            hypothesis <- gsub(tmp[w], newN[w], hypothesis, fixed=TRUE)
            cnames <- gsub(tmp[w], newN[w], cnames, fixed=TRUE)
        }
    cnames <- gsub(":",".", cnames)
    hypothesis <- gsub(":",".", hypothesis)
    ####
    
    l[cnames] <- 0
    chk <- grepl("=", hypothesis)
    R <- matrix(0, n,k)
    if (is.null(rhs))
    {
        hypothesis[!chk] <- paste(hypothesis[!chk],"=0",sep="")
    } else {
        if (any(chk))
            stop("hypothesis cannot contain = signs when rhs is not NULL")
        hypothesis <- paste(hypothesis, "=", rhs, sep="")
    }
    rhs <- numeric(n)
    tmp <- strsplit(hypothesis, "=")
    fl <- sapply(tmp, function(x) x[1])
    fr <- sapply(tmp, function(x) x[2])
    for (i in 1:n)
    {
        el <- parse(text=fl[i])
        er <- parse(text=fr[i])
        vr <- all.vars(er)
        vl <- all.vars(el)
        if (!all(c(vr,vl)%in%cnames))
            stop("wrong variable names. Special variable names may be a reason. Try to use a matrix R instead")
        if (length(vr) > 0)
        {
            tmp <- sapply(vr, function(v) try(eval(D(er, v), new.env()), silent=TRUE))
            if (!all(is.numeric(tmp)))
                stop("Bad hypothesis equations")
            R[i,match(vr,cnames)] <- -tmp
        }
        if (length(vl) > 0)
        {
            tmp <- sapply(vl, function(v) try(eval(D(el, v), new.env()), silent=TRUE))
            if (!all(is.numeric(tmp)))
                stop("Bad hypothesis equations")
            R[i,match(vl,cnames)] <- R[i,match(vl,cnames)] + tmp
        }
        rhs[i] <- eval(er, l) - eval(el, l)
    }
    list(R=R, rhs=rhs)
}

## print restriction on restricted models

setGeneric("printRestrict", function(object, ...)
    standardGeneric("printRestrict"))

setMethod("printRestrict", "rlinearModel",
          function(object){
              cst <- .printHypothesis(object@cstLHS, object@cstRHS, object@parNames)
              cat("Constraints:\n")
              for (i in 1:length(cst))
                  cat("\t", cst[i], "\n")
              cat("Restricted regression:\n\t")
              cat(.printRFct(object), "\n")
          })

setMethod("printRestrict", "rnonlinearModel",
          function(object){
              cat("Constraints:\n")
              for (i in 1:length(object@R))
                  {
                      cat("\t")
                      print(object@R[[i]])
                  }
          })

setMethod("printRestrict", "rformulaModel",
          function(object){
              cat("Constraints:\n")
              for (i in 1:length(object@R))
                  {
                      cat("\t")
                      print(object@R[[i]])
                  }
          })

setMethod("printRestrict", "rfunctionModel",
          function(object){
              cat("Constraints:\n")
              for (i in 1:length(object@R)) {
                  cat("\t")
                  print(object@R[[i]])
              }})

## print

setMethod("print", "rlinearModel",
          function(x)
          {
              callNextMethod()
              printRestrict(x)
          })

setMethod("print", "rformulaModel",
          function(x)
          {
              callNextMethod()
              printRestrict(x)
          })

setMethod("print", "rnonlinearModel",
          function(x)
          {
              callNextMethod()
              printRestrict(x)
          })

setMethod("print", "rfunctionModel",
          function(x) {
              callNextMethod()
              printRestrict(x)
          })

## restModel constructor

setGeneric("restModel", function(object, ...) standardGeneric("restModel"))

setMethod("restModel", signature("linearModel"),
          function(object, R, rhs=NULL)
          {
              if (is.character(R))
              {
                  res <- .makeHypothesis(object@parNames, R, rhs)
                  R <- res$R
                  rhs <- res$rhs
              } else {
                  if (is.null(rhs))
                      rhs <- rep(0,nrow(R))
              }
              res <- try(.imposeRestrict(R,rhs,object@parNames), silent=TRUE)
              if (any(class(res) == "try-error"))
                  stop("Failed to construct restricted model from the provided restrictions can you simplify it?")
              isEndo <- object@isEndo
              rtet <- res$theta
              res$isEndo <- c(crossprod(isEndo, rtet)) != 0
              new("rlinearModel",  cstLHS=R, cstRHS=rhs,
                  cstSpec=res, object)
          })

setMethod("restModel", signature("nonlinearModel"),
          function(object, R, rhs=NULL) {
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
              res <- .imposeNLRestrict(R, object)
              cstSpec <- list(newParNames = res$parNames,
                              originParNames=object@parNames,
                              k=res$k, theta0=res$theta0, fRHS=res$rhs, fLHS=res$lhs)
              new("rnonlinearModel", R=R, cstSpec=cstSpec, object)
          })

setMethod("restModel", signature("functionModel"),
          function(object, R, rhs=NULL) {
              if (!is.null(rhs))
                  warning("rhs is ignored for functional models")
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
              res <- .imposefRestrict(R, object)
              cstSpec <- list(newParNames = res$parNames,
                              originParNames=object@parNames,
                              k=res$k, theta0=res$theta0)
              new("rfunctionModel", R=R, cstSpec=cstSpec, object)
          })

setMethod("restModel", signature("formulaModel"),
          function(object, R, rhs=NULL) {
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
              res <- .imposeFORMRestrict(R, object)
              cstSpec <- list(newParNames = res$parNames,
                              originParNames=object@parNames,
                              k=res$k, theta0=res$theta0, fRHS=res$rhs, fLHS=res$lhs)
              new("rformulaModel", R=R, cstSpec=cstSpec, object)
          })

### Get the restriction matrices

setGeneric("getRestrict", function(object, ...)
    standardGeneric("getRestrict"))

setMethod("getRestrict", "rlinearModel",
          function(object, theta) {
              theta <- setCoef(as(object, "linearModel"), theta)
              R <- c(object@cstLHS%*%theta)
              cst <- .printHypothesis(object@cstLHS, object@cstRHS, object@parNames)
              list(dR=object@cstLHS, R=R, q=object@cstRHS, hypo=cst,
                   orig.R=object@cstLHS, orig.rhs=object@cstRHS)
          })

setMethod("getRestrict", "rnonlinearModel",
          function(object, theta) {
              dR <-numeric()
              R <- numeric()
              theta <- setCoef(as(object, "nonlinearModel"), theta)              
              for (r in object@R)
                  {
                      dlhs <- sapply(object@parNames, function(pn)
                          eval(D(r[[2]], pn), as.list(theta)))
                      drhs <- sapply(object@parNames, function(pn)
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

setMethod("getRestrict", "rformulaModel",
          function(object, theta) {
              getMethod("getRestrict", "rnonlinearModel")(object, theta)
          })


setMethod("getRestrict", "rfunctionModel",
          function(object, theta){
              getMethod("getRestrict", "rnonlinearModel")(object, theta)
          })

setMethod("getRestrict", "momentModel",
          function(object, theta, R, rhs=NULL) {
              robject <- restModel(object, R, rhs)
              getRestrict(robject, theta)
          })

## coef get the coefficients using the unrestricted representation

setMethod("coef", "rlinearModel",
          function(object, theta)
              {
                  cst <- object@cstSpec
                  if (length(theta)!=cst$k)
                      stop("Wrong number of coefficients")
                  if (cst$k == 0)
                      {
                          theta <- cst$minY
                          names(theta) <- object@parNames
                          return(theta)
                      }
                  tet <- cst$minY
                  tet2 <- apply(cst$theta, 1, function(x) sum(x*theta))
                  tet <- tet+tet2
                  names(tet) <- object@parNames
                  tet
              })

setMethod("coef", "rnonlinearModel",
          function(object, theta)
          {
              spec <- modelDims(object)
              theta <- setCoef(object, theta)
              theta2 <- rep(0,object@k)
              names(theta2) <- object@parNames
              theta2[names(theta)] <- theta
              chk <- sapply(object@R, function(r) is.numeric(r[[3]]))
              for (r in object@R[chk])
                  theta2[as.character(r[[2]])] <- r[[3]]
              for (r in object@R[!chk])
                  theta2[as.character(r[[2]])] <- eval(r[[3]], as.list(theta2))
              theta2
          })

setMethod("coef", "rfunctionModel",
          function(object, theta)
              getMethod("coef","rnonlinearModel")(object, theta)
          )

setMethod("coef", "rformulaModel",
          function(object, theta)
              getMethod("coef","rnonlinearModel")(object, theta)
          )


## Subsetting '['

setMethod("[", c("rfunctionModel", "numeric", "missing"),
          function(x, i, j){
              Call <- match.call(call=sys.call(sys.parent()))
              obj <- callNextMethod()
              obj@call <- Call
              obj
          })

## gmmfit

setMethod("gmmFit", signature("rlinearModel"), valueClass="gmmfit", 
          definition = function(model, type=c("twostep", "iter","cue", "onestep"),
              itertol=1e-7, initW=c("ident", "tsls"), weights="optimal", 
              itermaxit=100, efficientWeights=FALSE, ...) {
              Call <- try(match.call(call=sys.call(sys.parent())), silent=TRUE)
              if (inherits(Call,"try-error"))              
                  Call <- NULL
              cst <- model@cstSpec
              if (cst$k==0)
                  {
                      theta <- coef(model, numeric())
                      model <- as(model, "linearModel")                      
                      if (inherits(weights,"momentWeights"))
                          wObj <- weights
                      else
                          wObj <- evalWeights(model, theta=theta, w=weights)
                      obj <- evalGmm(model, theta, wObj)
                  } else {
                      obj <- callNextMethod()
                  }
              obj@call <- Call
              obj
          })

setMethod("gmmFit", signature("rnonlinearModel"), valueClass="gmmfit", 
          definition = function(model, type=c("twostep", "iter","cue", "onestep"),
              itertol=1e-7, initW=c("ident", "tsls"), weights="optimal", 
              itermaxit=100, efficientWeights=FALSE, theta0=NULL, ...) {
              Call <- try(match.call(call=sys.call(sys.parent())), silent=TRUE)
              if (inherits(Call,"try-error"))
                  Call <- NULL
              cst <- model@cstSpec
              if (cst$k==0)
                  {
                      theta <- coef(model, numeric())
                      model <- as(model, "nonlinearModel")                      
                      if (inherits(weights,"momentWeights"))
                          wObj <- weights
                      else
                          wObj <- evalWeights(model, theta=theta, w=weights)
                      obj <- evalGmm(model, theta, wObj, Call=FALSE)
                  } else {
                      obj <- callNextMethod()
                  }
              obj@call <- Call
              obj
          })

setMethod("gmmFit", signature("rformulaModel"), valueClass="gmmfit", 
          definition = function(model, type=c("twostep", "iter","cue", "onestep"),
              itertol=1e-7, initW=c("ident", "tsls"), weights="optimal", 
              itermaxit=100, efficientWeights=FALSE, theta0=NULL, ...) {
              Call <- try(match.call(call=sys.call(sys.parent())), silent=TRUE)
              if (inherits(Call,"try-error"))              
                  Call <- NULL
              cst <- model@cstSpec
              if (cst$k==0)
                  {
                      theta <- coef(model, numeric())
                      model <- as(model, "formulaModel")                      
                      if (inherits(weights,"momentWeights"))
                          wObj <- weights
                      else
                          wObj <- evalWeights(model, theta=theta, w=weights)
                      obj <- evalGmm(model, theta, wObj)
                  } else {
                      obj <- callNextMethod()
                  }
              obj@call <- Call
              obj              
          })


### momentStrength
### For now, there is no measure of moment strength in restricted models
### Have to figure out how to identify exluded instruments after
### the model has been modified.

setMethod("momentStrength", "rlinearModel",
          function(object, theta, vcovType = c("OLS", "HC", "HAC"))  {
              fstats <- NULL
              mess <- "No strength measure available for restricted models"  
              list(strength=fstats, mess=mess)
          })



### GEL specifics
######################

setMethod("gelFit", signature("rmomentModel"), valueClass="gelfit", 
          definition = function(model, gelType="EL", rhoFct=NULL,
                                initTheta=c("gmm", "modelTheta0"), theta0=NULL,
                                lambda0=NULL, vcov=FALSE, ...)
          {
              Call <- try(match.call(call=sys.call(sys.parent())), silent=TRUE)
              if (inherits(Call,"try-error"))
                  Call <- NULL
              k <- modelDims(model)$k
              if (k == 0)
                  return(evalGel(model, numeric(), ...))
              initTheta <- match.arg(initTheta)
              if (is.null(theta0))
              {
                  if (initTheta == "gmm")
                  {
                      theta0 <- gmmFit(model)@theta
                  } else {
                      theta0 <- modelDims(model)$theta0
                  }
              }
              obj <- getMethod("gelFit", "momentModel")(model=model, gelType=gelType,
                  rhoFct=rhoFct, initTheta=initTheta, theta0=theta0,
                  lambda0=lambda0, vcov=vcov, ...)
              obj@call <- Call
              obj
              })
