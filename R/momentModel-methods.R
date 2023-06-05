####### All methods with gmmModels (and its subclasses) signature
#################################################################

#######################  Print ########################
### The getGeneric for print is here only, so the file must be compiled
### before any other files containing print

setGeneric("print")

setMethod("print", "sSpec",
          function(x, digits=3, ...)
          {
              cat("Smoothing: ")
              if (x@kernel == "none")
              {
                  cat("none\n")
              } else {
                  cat(x@kernel, " kernel and ", sep = "")
                  cat(x@bwMet, " bandwidth", sep = "")
                  cat(" (", round(x@bw, digits), ")", sep = "")
              }
              cat("\n")
              invisible()
          })

setMethod("show", "sSpec", function(object) print(object))

setMethod("print", "momentModel",
          function(x, ...) {
              cat("Model based on moment conditions\n")
              cat("*********************************\n")
              cat("Moment type: ", strsplit(is(x)[1], "Model")[[1]][1], "\n", sep="")
              cat("Covariance matrix: ", x@vcov, sep="")
              if (x@vcov == "HAC")
                  {
                      cat(" with ", x@vcovOptions$kernel, " kernel and ", sep="")
                      if (is.numeric(x@vcovOptions$bw))
                          cat("Fixed  bandwidth (", round(x@vcovOptions$bw,3), ")",  sep="")
                      else
                          cat(x@vcovOptions$bw, " bandwidth",  sep="")
                  }
              if (x@smooth)
              {
                  cat("\n")
                  print(x@sSpec, ...)
              }
              if (x@vcov == "CL")
                  cat("\nClustered based on: ",
                      paste(colnames(x@vcovOptions$cluster), collapse=" and "), sep="")
              if (length(x@survOptions)>0)
                  cat("\nSurvey weights type: ", x@survOptions$type, sep="")
              cat("\n")
              d <- modelDims(x)
              cat("Number of regressors: ", d$k, "\n", sep="")
              cat("Number of moment conditions: ", d$q, "\n", sep="")
              if (!inherits(x, "functionGmm"))
                  cat("Number of Endogenous Variables: ", sum(x@isEndo), "\n", sep="")
              if (!is.null(x@survOptions$weights) && x@survOptions$type == "frequency")
                  cat("Implied sample size (sum of weights): ", d$n, "\n")
              else
                  cat("Sample size: ", d$n, "\n")
              invisible()
          })

setMethod("show", "momentModel", function(object) print(object))

### setCoef #######
### A generic for setting and validating the format of the coefficient

setGeneric("setCoef", function(model, ...) standardGeneric("setCoef"))

setMethod("setCoef", "momentModel",
          function(model, theta) {
              spec <- modelDims(model)
              if (length(theta) != length(spec$parNames))
                  stop(paste("Wrong number of coefficients (should be ",
                             spec$k, ")", sep=""))
              if (!is.null(names(theta)))
              {
                  chk <- !(names(theta)%in%spec$parNames)
                  if (any(chk))
                  {
                      mess <- paste("The following theta names are invalid: ",
                                    paste(names(theta)[chk], collapse=", ", sep=""),
                                    sep="")
                      stop(mess)
                  }
                  theta <- theta[match(spec$parNames, names(theta))]
              } else {
                  names(theta) <- spec$parNames
              }
              theta})

##### coef  ########
### For this, it only attach the names to theta

setMethod("coef", "momentModel",
          function(object, theta) {
              if (length(theta) != length(object@parNames))
                  stop("Wrong number of coefficients")
              if (!is.null(names(theta)))
              {
                  if (!all(names(theta)%in%object@parNames))
                      stop("theta has wrong names")
                  theta <- theta[match(object@parNames, names(theta))]
              } else {
                  names(theta) <- object@parNames
              }
              theta})

################## model.matrix and modelResponse #################
### I did not make model.response as generic because it is not
### a method in stats and I want different arguments

setGeneric("modelResponse", function(object, ...) standardGeneric("modelResponse"))

setMethod("modelResponse", signature("linearModel"),
          function(object)
          {
              model.response(object@modelF)
          })

setGeneric("model.matrix")
setMethod("model.matrix", signature("linearModel"),
          function(object, type=c("regressors","instruments"))
          {
              type <- match.arg(type)
              if (type == "regressors")
              {
                  ti <- attr(object@modelF, "terms")
                  mat <- as.matrix(model.matrix(ti, object@modelF)[,])
              } else {
                  ti <- attr(object@instF, "terms")
                  mat <- as.matrix(model.matrix(ti, object@instF)[,])
              }
              mat
          })

setMethod("model.matrix", signature("nonlinearModel"),
          function(object, type=c("regressors","instruments"))
          {
              type <- match.arg(type)
              if (type == "regressors")
              {
                  stop("no model.matrix of type regressors for nonlinear Gmm. set type to 'instruments' to get the matrix of instruments")
              } else {
                  ti <- attr(object@instF, "terms")
                  mat <- model.matrix(ti, object@instF)[,]
              }
              mat
          })

#################### residuals ###########################
### The getGeneric for residuals is here only, so the file must be compiled
### before any other files containing print


### rlinearGmm inherits from linearGmm
setGeneric("residuals")
setMethod("residuals", signature("linearModel"), function(object, theta){
    X <- model.matrix(object)
    Y <- modelResponse(object)
    e <- Y-c(X%*%theta)
    e
})

setMethod("residuals", signature("nonlinearModel"), 
          function(object, theta)
              {
                  res <- modelDims(object)
                  theta <- coef(object, setCoef(object, theta))
                  varList <- c(as.list(theta), as.list(object@modelF))
                  if (!is.null(res$fLHS))
                      {
                          lhs <- try(eval(res$fLHS, varList))
                          if (inherits(lhs, "try-error"))
                              stop("Cannot evaluate the LHS")
                      } else {
                          lhs <- 0
                      }
                  rhs <- try(eval(res$fRHS, varList))
                  if (inherits(rhs, "try-error"))
                      stop("Cannot evaluate the RHS")
                  c(lhs-rhs)
              })

################ evalMoment ##########################

setGeneric("evalMoment", function(object, theta, ...) standardGeneric("evalMoment"))

setMethod("evalMoment", signature("regModel"),
          function(object, theta) {
              e <- residuals(object, theta)
              Z <- model.matrix(object, "instruments")
              gt <- as.matrix(Z)*e
              if (object@smooth)
                  gt <- stats::kernapply(gt, object@sSpec@w)
              gt
          })

setMethod("evalMoment", signature("functionModel"),
          function(object, theta) {
              theta <- coef(object, setCoef(object, theta))              
              gt <- object@fct(theta, object@X)
              if (!is.null(sub <- attr(object@X, "subset")))
                  gt <- gt[,sub]
              if (object@smooth)
                  gt <- stats::kernapply(gt, object@sSpec@w)
              gt
          })

setMethod("evalMoment", signature("formulaModel"),
          function(object, theta) {
              res <- modelDims(object)
              theta <- coef(object, setCoef(object, theta))              
              varList <- c(as.list(theta), as.list(object@modelF))
              gt <- sapply(1:res$q, function(i) {
                  if (!is.null(res$fLHS[[i]]))
                  {
                      lhs <- try(eval(res$fLHS[[i]], varList))
                      if (inherits(lhs, "try-error"))
                          stop("Cannot evaluate the LHS")
                  } else {
                      lhs <- 0
                  }
                  if (!is.null(res$fRHS[[i]]))
                  {
                      rhs <- try(eval(res$fRHS[[i]], varList))
                      if (inherits(rhs, "try-error"))
                          stop("Cannot evaluate the RHS")
                  } else {
                      lhs <- 0
                  }
                  c(lhs-rhs)})
              if (object@smooth)
                  gt <- stats::kernapply(gt, object@sSpec@w)
              gt
          })


################ evalDresiduals ##########################

setGeneric("Dresiduals", function(object, theta, ...) standardGeneric("Dresiduals"))

setMethod("Dresiduals", signature("linearModel"),
          function(object, theta) {
              -model.matrix(object)
          })

setMethod("Dresiduals", signature("nonlinearModel"),
          function(object, theta) {
              theta <- coef(object, setCoef(object, theta))              
              res <- modelDims(object)
              varList <- c(as.list(theta), as.list(object@modelF))
              De <- numeric()
              for (i in names(theta))
                  {
                      if (!is.null(res$fLHS))
                          d <- eval(D(res$fLHS, i), varList)      
                      else
                          d <- 0
                      De <-  cbind(De, d-matrix(eval(D(res$fRHS, i), varList),res$n,1))
                  }
              colnames(De) <- object@parNames
              De
          })

############### modelDims #######################

setGeneric("modelDims", function(object, ...) standardGeneric("modelDims"))

setMethod("modelDims", "linearModel",
          function(object) {
              n <- if (object@smooth)
                       object@n-2*object@sSpec@w$m
                   else
                       object@n
              list(k=object@k, q=object@q, n=n, parNames=object@parNames,
                   momNames=object@momNames, isEndo=object@isEndo)
          })

setMethod("modelDims", "nonlinearModel",
          function(object) {
              n <- if (object@smooth)
                       object@n-2*object@sSpec@w$m
                   else
                       object@n              
              list(k=object@k, q=object@q, n=n, parNames=object@parNames,
                   momNames=object@momNames, theta0=object@theta0,
                   fRHS=object@fRHS, fLHS=object@fLHS, isEndo=object@isEndo)
          })

setMethod("modelDims", "functionModel",
          function(object) {
              n <- if (object@smooth)
                       object@n-2*object@sSpec@w$m
                   else
                       object@n              
              list(k=object@k, q=object@q, n=n, parNames=object@parNames,
                   momNames=object@momNames, theta0=object@theta0,
                   fct=object@fct, dfct=object@dfct, isEndo=object@isEndo)
          })

setMethod("modelDims", "formulaModel",
          function(object) {
              n <- if (object@smooth)
                       object@n-2*object@sSpec@w$m
                   else
                       object@n              
              list(k=object@k, q=object@q, n=n, parNames=object@parNames,
                   momNames=object@momNames, theta0=object@theta0,
                   fRHS=object@fRHS, fLHS=object@fLHS, isEndo=object@isEndo,
                   isMDE=object@isMDE)
          })

################ evalDMoment ##########################

setGeneric("evalDMoment", function(object, ...) standardGeneric("evalDMoment"))

setMethod("evalDMoment", signature("regModel"),
          function(object, theta, impProb=NULL, lambda=NULL)
          {
              De <- Dresiduals(object, theta)
              Z <- model.matrix(object, "instrument")
              spec <- modelDims(object)
              if (is.null(impProb))
                  impProb <- 1/spec$n
              if (!is.null(lambda))
              {
                  if (length(lambda) != spec$q)
                      stop("The length of lambda must be equal to the number of conditions")
                  Z <- De*c(Z%*%lambda)
                  if (object@smooth)
                      Z <- stats::kernapply(Z, object@sSpec@w)
                  G <- Z*impProb
                  colnames(G) <- spec$parNames
                  rownames(G) <- NULL
              } else {
                  G <- apply(De,2, function(x)
                  {
                      tmp <- Z*x
                      if (object@smooth)
                          tmp <- stats::kernapply(tmp, object@sSpec@w)
                      colSums(tmp*impProb)
                  })
                  if (!is.matrix(G))
                      G <- matrix(G,  spec$q, spec$k)
                  dimnames(G) <- list(spec$momNames, spec$parNames)
              }
              G
          })

setMethod("evalDMoment", signature("functionModel"),
          function(object, theta, impProb=NULL, lambda=NULL)
          {
              spec <- modelDims(object)
              theta <- coef(object, setCoef(object, theta))              
              if (object@smooth && !is.null(object@dfct))
              {
                  object@dfct <- NULL
                  warning("Cannot provide dfct for smoothed models. dfct set to NULL")
              }              
              if (!is.null(lambda))
              {
                  if (length(lambda) != spec$q)
                      stop("The length of lambda must be equal to the number of conditions")
                  if (is.null(impProb))
                      impProb <- 1/spec$n              
                  f_lam <- function(theta, object, lambda)
                  {
                      gt <- evalMoment(object, theta)
                      c(gt%*%lambda)
                  }
                  env <- new.env()
                  assign("theta", theta, envir=env)
                  assign("object", object, envir=env)
                  assign("lambda", lambda, envir=env)
                  assign("f_lam", f_lam, envir=env)
                  G <- numericDeriv(quote(f_lam(theta, object, lambda)), "theta", env)
                  G <- attr(G, "gradient")*c(impProb)
                  colnames(G) <- spec$parNames
                  rownames(G) <- NULL
                  return(G)
              } 
              if (is.null(object@dfct))
              {
                  f <- function(theta, object, impProb)
                  {
                      gt <- evalMoment(object, theta)
                      if (is.null(impProb))
                          colMeans(gt)
                      else
                          colSums(gt*impProb)
                  }
                  env <- new.env()
                  assign("theta", theta, envir=env)
                  assign("object", object, envir=env)
                  assign("impProb", impProb, envir=env)
                  assign("f", f, envir=env)
                  G <- numericDeriv(quote(f(theta, object, impProb)), "theta", env)
                  G <- attr(G, "gradient")
              } else {
                  G <- object@dfct(theta, object@X)
              }
              if (!is.matrix(G))
                  G <- matrix(G,  spec$q, spec$k)
              dimnames(G) <- list(spec$momNames, spec$parNames)
              G
          })

setMethod("evalDMoment", signature("formulaModel"),
          function(object, theta, impProb=NULL, lambda=NULL)
          {
              theta <- coef(object, setCoef(object, theta))              
              spec <- modelDims(object)              
              if (is.null(impProb))
                  impProb <- 1/spec$n
              varList <- c(as.list(theta), as.list(object@modelF))
              if (!is.null(lambda))
              {
                  if (length(lambda) != spec$q)
                      stop("The length of lambda must be equal to the number of conditions")
                  f <- function(theta, lambda, eq)
                  {
                      Gi <- sapply(1:spec$q, function(j)
                      {
                          if (!is.null(eq[[j]]))
                          {
                              tmp <- eval(D(eq[[j]], i), varList)*lambda[j]
                              if (length(tmp)==1)
                                  tmp <- rep(tmp, spec$n)
                              if (object@smooth)
                                  tmp <- stats::kernapply(tmp, object@sSpec@w)
                              tmp <- c(tmp*impProb)
                          } else {
                              tmp <- numeric(spec$n)
                          }
                          tmp
                      })
                      rowSums(Gi)
                  }
                  nG <- list(NULL, spec$parNames)
              } else {
                  f <- function(theta, lambda, eq)
                  {
                      Gi <- sapply(1:spec$q, function(j)
                      {                          
                          if (!is.null(eq[[j]]))
                          {
                              tmp <- eval(D(eq[[j]], i), varList)
                              if (length(tmp)>1)
                              {
                                  if (object@smooth)
                                      tmp <- stats::kernapply(tmp, object@sSpec@w)
                                  tmp <- sum(tmp*impProb)
                              }
                          } else {
                              tmp <- 0
                          }
                          c(tmp)
                      })
                      Gi
                  }
                  nG <- list(spec$momNames, spec$parNames)
              }                  
              G <- numeric()
              for (i in names(theta))
              {
                  lhs <- f(i, lambda, spec$fLHS)
                  rhs <- f(i, lambda, spec$fRHS)
                  G <- cbind(G, lhs-rhs)
              }
              if (!is.matrix(G))
                  G <- matrix(G,  spec$q, spec$k)
              dimnames(G) <- nG 
              G
          })


###########   estfun :  Don't like it ###############

estfun.momentFct <- function(x, ...) x


##########  vcovHAC from sandwich #################

setMethod("vcovHAC", "momentModel",
          function (x, theta) { 
              if (x@vcov != "HAC")
              {
                  warning("Model set as ", x@vcov, ". The default HAC options are used")
                  x@vcov <- "HAC"
                  x@sSpec <- new("sSpec")
                  x@smooth <- FALSE
                  x@vcovOptions <- .getVcovOptions("HAC")
              }
              gmat <- evalMoment(x, theta)              
              if (x@centeredVcov) 
                  gmat <- scale(gmat, scale = FALSE)
              class(gmat) <- "momentFct"
              options <- x@vcovOptions
              if (is.character(options$bw))
              {
                  if (options$bw == "Wilhelm")
                  {
                      G <- evalDMoment(x, theta)
                      obj <- list(gt = gmat, G = G)
                      class(obj) <- "gmm"
                  } else {
                      obj <- gmat
                  }
                  bwFct  <- get(paste("bw",options$bw,sep=""))
                  bwArgs <- options
                  bwArgs$bw <- NULL
                  bwArgs$tol <- NULL
                  bwArgs$x <- obj
                  bw <- do.call(bwFct, bwArgs)
              } else {
                  bw <- options$bw
              }
              weights <- weightsAndrews(x = gmat, bw = bw, kernel = options$kernel, 
                                        prewhite = options$prewhite, tol = options$tol)
              w <- vcovHAC(x = gmat, order.by = NULL, weights = weights, 
                           prewhite = options$prewhite, sandwich = FALSE,
                           ar.method = options$ar.method)
              attr(w, "Spec") <- list(weights = weights, bw = bw, kernel = options$kernel)
              w
          })

################ vcov  ##########################

setGeneric("vcov")

setMethod("vcov", signature("momentModel"),
          function(object, theta){              
              if ((inherits(object, "functionModel") || inherits(object, "formulaModel")) &
                  object@vcov == "iid")
                  object@vcov <- "MDS"
              n <- modelDims(object)$n
              if (object@vcov == "MDS")
                  {
                      gt <- evalMoment(object, theta)
                      if (object@centeredVcov)
                          gt <- scale(gt, scale=FALSE)
                      w <- crossprod(gt)/n
                      if (object@smooth)
                          w <- w*object@sSpec@bw/object@sSpec@k[2]
                  } else if (object@vcov == "iid") {
                      sig <- sd(residuals(object, theta))
                      Z <- model.matrix(object, "instrument")
                      w <- sig^2*crossprod(Z)/nrow(Z)
                  } else if (object@vcov == "CL") {
                      gt <- evalMoment(object, theta)
                      class(gt) <- "momentFct"
                      opt <- object@vcovOptions
                      opt$x <- gt
                      w <- do.call(meatCL, opt)                                            
                  } else {
                      w <- vcovHAC(object, theta)
                  }
              w})

################### weights Object and methods: Is it too much??? #################


setGeneric("evalWeights", function(object, ...) standardGeneric("evalWeights"))

setMethod("evalWeights", signature("momentModel"),valueClass="momentWeights",
          function(object, theta=NULL, w="optimal", ...)
          {
              wSpec <- list()
              if (is.matrix(w))
              {
                  type <- "weights"
              } else {
                  if (w == "ident")
                  {
                      type <- "weights"
                  } else {
                      if (inherits(object, c("formulaModel", "functionModel")) &
                          object@vcov == "iid")
                          object@vcov <- "MDS"
                      if (object@vcov == "MDS")
                      {
                          gt <- evalMoment(object, theta)
                          if (object@centeredVcov)
                              gt <- scale(gt, scale=FALSE)
                          w <- qr(gt/sqrt(nrow(gt)))
                          if (w$rank < object@q)
                              warning("The moment matrix is not full column rank")
                          type <- "qr"
                      } else if (object@vcov == "iid") {
                          sig <- mean(residuals(object, theta)^2)
                          sig <- sqrt(sig)
                          ti <- attr(object@instF, "terms")
                          Z <- model.matrix(ti, object@instF)[,]
                          w <- qr(sig*Z/sqrt(nrow(Z)))
                          if (w$rank < object@q)
                              warning("The moment matrix is not full column rank")
                          type <- "qr"
                      } else if (object@vcov == "CL") {
                          gt <- evalMoment(object, theta)
                          class(gt) <- "momentFct"
                          opt <- object@vcovOptions
                          opt$x <- gt
                          w <- chol(do.call(meatCL, opt))
                          type <- "chol"
                      } else {
                          w <- vcovHAC(object, theta)
                          wSpec <- attr(w,"Spec")
                          w <- chol(w[,])
                          type <- "chol"
                      }
                  }
              }
              new("momentWeights", type=type, w=w, wSpec=wSpec)
          })

############ evalGmmObj #################################

setGeneric("evalGmmObj", function(object, theta, wObj, ...)
    standardGeneric("evalGmmObj"))

setMethod("evalGmmObj", signature("momentModel", "numeric", "momentWeights"),
          function(object, theta, wObj, ...)
              {
                  gt <- evalMoment(object, theta)
                  n <- modelDims(object)$n
                  gt <- colMeans(gt)
                  obj <- quadra(wObj, gt)
                  n*obj
              })

#########################  solveGmm  #########################

setGeneric("solveGmm", function(object, wObj, ...) standardGeneric("solveGmm"))

setMethod("solveGmm", signature("linearModel", "momentWeights"),
          function(object, wObj, theta0=NULL, ...)
          {
              X <- model.matrix(object)
              Z <- model.matrix(object, "instrument")
              Y <- modelResponse(object)
              d <- modelDims(object)
              n <- d$n
              Sig.zy <- crossprod(Z,Y)/n
              Sig.zx <- crossprod(Z,X)/n
              if (d$q == d$k)
              {
                  T1 <- Sig.zx
                  T2 <- Sig.zy
              } else {
                  T1 <- quadra(wObj, Sig.zx)
                  T2 <- quadra(wObj, Sig.zx, Sig.zy)
              }
              theta <- c(solve(T1, T2))
              names(theta) <- d$parNames
              list(theta=theta, convergence=NULL)
          })

setMethod("solveGmm", signature("allNLModel", "momentWeights"),
          function(object, wObj, theta0=NULL, algo=c("optim","nlminb"), ...)
          {
                  algo <- match.arg(algo)
                  if (is.null(theta0))
                      theta0 <- modelDims(object)$theta0
                  g <- function(theta, wObj, object)
                      evalGmmObj(object, theta, wObj)
                  dg <- function(theta, wObj, object)
                      {
                          gt <- evalMoment(object, theta)
                          n <- nrow(gt)
                          gt <- colMeans(gt)
                          G <- evalDMoment(object, theta)
                          obj <- 2*n*quadra(wObj, G, gt)
                          obj
                      }
                  if (algo == "optim")
                      {
                          if ("method" %in% names(list(...)))
                              res <- optim(par=theta0, fn=g, gr=dg, 
                                           object=object, wObj=wObj, ...)
                          else
                              res <- optim(par=theta0, fn=g, gr=dg, method="BFGS",
                                           object=object, wObj=wObj, ...)
                      } else {
                          res <- nlminb(start=theta0, objective=g, gradient=dg,
                                        object=object, wObj=wObj, ...)
                      }
                  theta <- res$par
                  names(theta) <- modelDims(object)$parNames
                  list(theta=theta, convergence=res$convergence)
              })


##################### momentStrength ####################

setGeneric("momentStrength", function(object, ...) standardGeneric("momentStrength"))

setMethod("momentStrength", signature("nonlinearModel"), 
          function(object, theta=NULL, ...)
          {
              list(strength=NULL, mess=NULL)
          })

setMethod("momentStrength", signature("functionModel"), 
          function(object, theta=NULL, ...)
          {
              list(strength=NULL, mess=NULL)
          })

setMethod("momentStrength", signature("formulaModel"), 
          function(object, theta=NULL, ...)
          {
              list(strength=NULL, mess=NULL)
          })

setMethod("momentStrength", signature("linearModel"), 
          function(object, theta, vcovType=c("OLS","HC","HAC","CL")){
              spec <- modelDims(object)
              getF <- function(i)
              {
                  resu <- lm(X[, i] ~ Z - 1)
                  v <- switch(vcovType, OLS = vcov(resu),
                              HC = vcovHC(resu, "HC1"),
                              HAC = vcovHAC(resu),
                              CL = do.call(vcovCL, c(object@vcovOptions, list(x = resu))))
                  v <- v[!exoInst, !exoInst]
                  b <- coef(resu)[!exoInst]
                  f <- b %*% solve(v, b)/df1
                  df2 <- resu$df.residual
                  c(f, df1, df2)
              }
              EndoVars <- !(spec$parNames %in% spec$momNames)
              exoInst <- spec$momNames %in% spec$parNames
              vcovType <- match.arg(vcovType)
              if (all(!EndoVars)) {
                  fstats <- NULL
                  mess <- "No endogenous variables: no strength measure"
              }
              else {
                  X <- model.matrix(object)
                  X <- X[, EndoVars, drop = FALSE]
                  Z <- model.matrix(object, "instrument")
                  fstats <- matrix(ncol = 0, nrow = 3)
                  df1 <- sum(!exoInst)
                  mess <- "Instrument strength based on the F-Statistics of the first stage OLS"
                  addM <- FALSE
                  for (i in 1:ncol(X)) {
                      tmp <- try(getF(i), silent=TRUE)
                      if (inherits(tmp, "try-error"))
                      {
                          fstats <- cbind(fstats, c(NA, NA, NA))
                          addM <- TRUE
                      } else {
                          fstats <- cbind(fstats, tmp)
                      }
                  }
                  if (addM)
                      mess <- paste(mess, "\n",
                                    "(Failed to compute some first stage F-statistics)",
                                    sep="")
                  fstats <- rbind(fstats, 1 - pf(fstats[1, ], fstats[2,], fstats[3, ]))
                  colnames(fstats) <- colnames(X)
                  rownames(fstats) <- c("Stats", "df1", "df2", "pv")
                  fstats <- t(fstats)

              }
              list(strength = fstats, mess = mess)
          })

### Subsetting models

setGeneric("[")
setMethod("[", c("regModel", "numeric", "missing"),
          function(x, i, j){
              i <- unique(as.integer(i))
              spec <- modelDims(x)
              q <- spec$q
              if (!all(abs(i) %in% (1:q))) 
                  stop("SubMoment must be between 1 and q")
              momNames <- x@momNames[i]
              if (length(momNames)<spec$k)
                  stop("The model is under-identified")
              if (momNames[1] == "(Intercept)") 
                  f <- reformulate(momNames[-1], NULL, TRUE)
              else 
                  f <- reformulate(momNames, NULL, FALSE)
              instF <- model.frame(f, x@instF)
              x@q <- length(momNames)
              x@instF <- instF
              x@momNames <- momNames
              x
          })


setMethod("[", c("functionModel", "numeric", "missing"),
          function(x, i, j){
              i <- unique(as.integer(i))
              spec <- modelDims(x)
              q <- spec$q
              if (!all(abs(i) %in% (1:q))) 
                  stop("SubMoment must be between 1 and q")
              if (length(i)==q)
                  return(x)
              momNames <- x@momNames[i]
              if (length(momNames)<spec$k)
                  stop("The model is under-identified")
              attr(x@X, "subset") <- i
              x@q <- length(momNames)
              x@momNames <- momNames
              x
          })

setMethod("[", c("formulaModel", "numeric", "missing"),
          function(x, i, j){
              i <- unique(as.integer(i))
              spec <- modelDims(x)
              q <- spec$q
              if (!all(abs(i) %in% (1:q))) 
                  stop("SubMoment must be between 1 and q")
               if (length(i)==q)
                  return(x)
               momNames <- x@momNames[i]
               if (length(momNames)<spec$k)
                   stop("The model is under-identified")
              x@fRHS <- x@fRHS[i]
              x@fLHS <- x@fLHS[i]
              x@q <- length(momNames)
              x@momNames <- momNames
              x
           })

setMethod("[", c("momentModel", "missing", "missing"),
          function(x, i, j) x)

### Observation subset

setGeneric("subset")
setMethod("subset", "regModel",
          function(x, i) {
              x@modelF <- x@modelF[i,,drop=FALSE]
              x@instF <- x@instF[i,,drop=FALSE]
              if (!is.null(x@vcovOptions$cluster))
                  x@vcovOptions$cluster <- x@vcovOptions$cluster[i,,drop=FALSE]
              if (!is.null(x@survOptions$weights))
                  x@survOptions$weights <- x@survOptions$weights[i]
              x@n <- nrow(x@modelF)
              x})

setMethod("subset", "functionModel",
          function(x, i) {
              if (is.matrix(x@X) || is.data.frame(x@X))
                  x@X <- x@X[i,,drop=FALSE]
              else if (is.numeric(x@X))
                  x@X <- x@X[i, drop=FALSE]
              else
                  stop("X is not subsetable")
              if (!is.null(x@vcovOptions$cluster))
                  x@vcovOptions$cluster <- x@vcovOptions$cluster[i,,drop=FALSE]
              if (!is.null(x@survOptions$weights))
                  x@survOptions$weights <- x@survOptions$weights[i]              
              x@n <- NROW(x@X)
              x})

setMethod("subset", "formulaModel",
          function(x, i) {
              x@modelF <- x@modelF[i,,drop=FALSE]
              if (!is.null(x@vcovOptions$cluster))
                  x@vcovOptions$cluster <- x@vcovOptions$cluster[i,,drop=FALSE]
              if (!is.null(x@survOptions$weights))
                  x@survOptions$weights <- x@survOptions$weights[i]              
              x@n <- nrow(x@modelF)
              x})

## modelFits

setGeneric("gmmFit", function(model, ...) standardGeneric("gmmFit"))

setMethod("gmmFit", signature("formulaModel"), valueClass="gmmfit", 
          definition = function(model, type=c("twostep", "iter","cue", "onestep"),
                                itertol=1e-7, initW=c("ident", "tsls"), weights="optimal", 
                                itermaxit=100, efficientWeights=FALSE, theta0=NULL, ...)
          {
              Call <- try(match.call(call=sys.call(sys.parent())), silent=TRUE)
              if (inherits(Call,"try-error"))
                  Call <- NULL
              if (model@isMDE && model@centeredVcov)
              {
                  if (is.character(weights) && weights == "optimal")
                  {
                      spec <- modelDims(model)
                      wObj <- evalWeights(model, spec$theta0, "optimal")
                      met <- getMethod("gmmFit", "momentModel")
                      res <- met(model, weights=wObj, efficientWeights=TRUE, ...)
                      res@type <- "mde"
                  } else {
                      res <- callNextMethod()
                  }
              } else {
                  res <- callNextMethod()
              }
              res@call <- Call              
              return(res)
          })

setMethod("gmmFit", signature("momentModel"), valueClass="gmmfit", 
         definition = function(model, type=c("twostep", "iter","cue", "onestep"),
              itertol=1e-7, initW=c("ident", "tsls"), weights="optimal", 
              itermaxit=100, efficientWeights=FALSE, theta0=NULL, ...)
         {
             Call <- try(match.call(call=sys.call(sys.parent())), silent=TRUE)
             if (inherits(Call,"try-error"))
                 Call <- NULL
             chk <- validObject(model)                  
             type <- match.arg(type)
             initW <- match.arg(initW)
             i <- 1L
             chk <- validObject(model, TRUE)
             if (!chk)
                 stop("model is not a valid momentModel object")
             if (initW == "tsls" && !is(model,"linearModel"))
                 stop("initW='tsls' is for linear models only")
             if (is.character(weights) && !(weights%in%c("optimal","ident")))
                 stop("weights is a matrix or one of 'optimal' or 'ident'")
             spec <- modelDims(model)
             if (spec$q==spec$k)
             {
                 ## This allow to weight the moments in case of
                 ## large scale difference.
                 if (!is.matrix(weights) && !inherits(weights,"momentWeights"))
                     weights <- "ident"
                 type <- "onestep"
             } else if (type == "onestep" && !is.matrix(weights)) {
                 weights <- "ident"
             } else if (is.matrix(weights) || inherits(weights,"momentWeights")) {
                 type <- "onestep"
             } else if (weights == "ident") {
                 type <- "onestep"
             }
             if (type == "onestep")
             {
                 if (inherits(weights,"momentWeights"))
                     wObj <- weights
                 else
                     wObj <- evalWeights(model, w=weights)
                 res <- solveGmm(model, wObj, theta0, ...)
                 convergence <- res$convergence
                 efficientGmm <- ifelse(is.character(weights), FALSE,
                                        efficientWeights)
                 ans <- new("gmmfit", theta=res$theta,
                            convergence=convergence, convIter=NULL, type=type,
                            wObj=wObj, model=model, call=Call, niter=i,
                            efficientGmm=efficientGmm)
                 return(ans)
             }
             if (is(model,"linearModel"))
             {
                 if (model@vcov == "iid")
                     if (is.character(weights) && weights == "optimal")
                     {
                         res <- tsls(model)
                         res@call <- Call
                         return(res)
                     }
             }
             if (type == "twostep")
             {
                 itermaxit <- 1
             }
             if (initW=="tsls")
             {                          
                 theta0 <- coef(tsls(model))
             } else {
                 wObj <- evalWeights(model, NULL, "ident")
                 theta0 <- solveGmm(model, wObj, theta0, ...)$theta
             }
             bw <- model@vcovOptions$bw
             if (type != "cue")
             {
                 while(TRUE)
             {
                 wObj <- evalWeights(model, theta0, "optimal")
                 if (model@vcov=="HAC" && is.character(bw))
                     model@vcovOptions$bw <- wObj@wSpec$bw
                 res <- solveGmm(model, wObj, theta0, ...)
                 theta1 <- res$theta
                 convergence <- res$convergence
                 crit <- sqrt( sum((theta1-theta0)^2)/(1+sqrt(sum(theta0^2))))
                 if (crit < itertol & type=="iter")
                 {
                     convIter <- 0
                     break
                 }
                 i <- i + 1L
                 theta0 <- theta1
                 if (i>itermaxit)
                 {
                     if (type=="twostep")
                         convIter <- NULL
                     else
                         convIter <- 1
                     break                                      
                 }                              
             }      
             } else {
                 convIter <- NULL
                 if (model@vcov=="HAC" && is.character(bw))
                 {
                     w <- vcov(model, theta0)
                     model@vcovOptions$bw <- attr(w, "Spec")$bw
                 }
                 obj <- function(theta, model)
                 {
                     wObj <- evalWeights(model, theta, "optimal")
                     evalGmmObj(model, theta, wObj)
                 }
                 res <- optim(theta0, obj, model=model,
                              ...)
                 theta1 <- res$par
                 convergence <- res$convergence
                 wObj <- evalWeights(model, theta1, "optimal")                 
             }
             model@vcovOptions$bw <- bw
             names(theta1) <- spec$parNames
             new("gmmfit", theta=theta1, convergence=convergence, type=type,
                 wObj=wObj, model=model, convIter=convIter, call=Call,
                 niter=i, efficientGmm=TRUE)
         })

## tsls

setGeneric("tsls", function(model, ...) standardGeneric("tsls"))

setMethod("tsls", signature("linearModel"), valueClass="tsls", 
          function(model)
          {
              Call <- try(match.call(call=sys.call(sys.parent())), silent=TRUE)
              if (inherits(Call,"try-errors"))
                  Call <- NULL
              chk <- validObject(model)
              X <- model.matrix(model)
              Z <- model.matrix(model, "instrument")
              Y <- modelResponse(model)
              spec <- modelDims(model)
              EndoVars <- !(spec$parNames %in% spec$momNames)
              if (any(EndoVars))
                  {
                      res <- lm(X[,EndoVars]~Z)
                      X[,EndoVars] <- fitted(res)
                  }
              theta <- lm.fit(X, Y)$coefficients
              names(theta) <- spec$parNames
              vcov <- model@vcov
              model@vcov <- "iid"
              efficientGmm <- vcov == "iid"
              wObj <- evalWeights(model, theta, "optimal")
              model@vcov <- vcov
              obj <- new("tsls", theta=theta, convergence=NULL, type="tsls",
                         wObj=wObj, model=model, convIter=NULL, call=Call,
                         niter=1L, efficientGmm=efficientGmm)
              obj
          })


## evalModels

setGeneric("evalGmm", function(model, ...) standardGeneric("evalGmm"))

setMethod("evalGmm", signature("momentModel"),
          function(model, theta, wObj=NULL, ...) {
              spec <- modelDims(model)
              Call <- try(match.call(call=sys.call(sys.parent())), silent=TRUE)
              if (inherits(Call,"try-error"))
                  Call <- NULL
              theta <- setCoef(model, theta)
              if (is.null(wObj))
                  wObj <- evalWeights(model, theta)
              new("gmmfit", theta=theta, convergence=NULL, convIter=NULL,
                  call=Call, type="eval", wObj=wObj, niter=0L, efficientGmm=FALSE,
                  model=model)
          })

## update

setGeneric("update")

setMethod("update", "list",
          function(object, ...)
          {              
              if (length(object) == 0)
                  return(object)
              if (any(names(object) == ""))
                  stop("cannot update a list if some elements are nameless")
              arg <- list(...)
              arg <- arg[na.omit(match(names(object), names(arg)))]
              for (n in names(arg))
                  object[[n]] <- arg[[n]]
              object
          })

setMethod("update", "momentModel",
          function(object, ...)
          {
              arg <- list(...)
              allowed <- c("vcov","vcovOptions", "centeredVcov",
                           "survOptions", "smooth")              
              arg <- arg[na.omit(match(allowed, names(arg)))]              
              if (length(arg) == 0)                  
                  return(object)
              chk <- FALSE
              if (!is.null(arg$smooth))
              {
                  if (!is.logical(arg$smooth))
                      stop("smooth must be logical")
                  if (length(arg$smooth) != 1L)
                      stop("smooth cannot be a vector")
                  if (object@smooth != arg$smooth)
                  {
                      if (!arg$smooth)
                      {
                          object@smooth <- FALSE
                          object@sSpec <- new("sSpec")
                      } else {
                          object@vcov <- "MDS"
                          object@smooth <- TRUE
                          chk <- TRUE
                      }
                  }
              }
              if (!is.null(arg[["vcov"]]) && !object@smooth)
                  object@vcov <- arg[["vcov"]]
              if (object@vcov == "HAC" || object@smooth)
              {                  
                  if (is.null(arg$vcovOptions))
                      arg$vcovOptions <- list()
                  if (length(object@vcovOptions))
                  {
                      tmp <- c(arg$vcovOptions, list(object=object@vcovOptions))
                      arg$vcovOptions <- do.call(update, tmp)
                  }
                  arg$vcovOptions <- .getVcovOptions(object@vcov, NULL, arg$vcovOptions,
                                                     object@smooth)              
                  if (object@smooth && !identical(arg$vcovOptions, object@vcovOptions))
                      chk <- TRUE
                  object@vcovOptions <- arg$vcovOptions
              }
              if (is.null(arg$survOptions))
                  arg$survOptions <- list()
              object@survOptions <- .getSurvOptions(NULL, arg$survOptions)              
              if (!is.null(arg$centeredVcov))
                  object@centeredVcov <- arg$centeredVcov
              if (chk)
                  object@sSpec <- kernapply(object, smooth=FALSE)
              object
              })

## kernapply
        
setGeneric("kernapply")

setMethod("kernapply", "momentModel",
          function(x, theta=NULL, smooth=TRUE, ...)
          {
              sSpec <- x@sSpec
              x@smooth <- FALSE
              x@sSpec <- new("sSpec")              
              if (smooth) {
                  if (is.null(theta))
                      stop("theta0 is needed to compute the smoothed moments")
                  gt <- evalMoment(x, theta)
                  sx <- stats::kernapply(gt, sSpec@w)
                  ans <- list(smoothx = sx, w = sSpec@w,
                              bw = sSpec@bw, k = sSpec@k)
                  return(ans)
              }
              if (is.null(theta))        
                  theta <- gmmFit(x, weights="ident")@theta
              gt <- evalMoment(x, theta)
              gt <- scale(gt, scale=FALSE)
              class(gt) <- "momentFct"
              vspec <- x@vcovOptions
              if (!(vspec$kernel%in%c("Bartlett","Parzen")))
                  x@vcovOptions$kernel <- "Bartlett"
              kernel <- switch(x@vcovOptions$kernel,
                               Bartlett="Truncated",
                               Parzen="Bartlett")
              k <- switch(kernel,
                          Truncated=c(2,2),
                          Bartlett=c(1,2/3))
              if (is.character(vspec$bw))
              {
                  bw <- get(paste("bw", vspec$bw, sep = ""))
                  bw <- bw(gt, kernel = vspec$kernel, prewhite = vspec$prewhite,
                           ar.method = vspec$ar.method, approx = vspec$approx)
                  bwMet <- vspec$bw
              } else {
                  bw <- x@vcovOptions$bw
                  bwMet <- "Fixed"
              } 
              w <- weightsAndrews(gt, bw = bw, kernel = kernel, prewhite = vspec$prewhite, 
                                  ar.method = vspec$ar.method, tol = vspec$tol,
                                  verbose = FALSE, approx = vspec$approx)
              rt <- length(w)
              if (rt >= 2)
              {
                  w <- c(w[rt:2], w)
                  w <- w/sum(w)
                  w <- kernel(w[rt:length(w)])
              } else {
                  w <- kernel(1)
              }
              new("sSpec", k=k, kernel=kernel, bw=bw, w=w, bwMet=bwMet)
          })

###################  GEL specifics ################
###################################################

############ evalGelObj #################################

setGeneric("evalGelObj", function(object, theta, lambda, ...)
    standardGeneric("evalGelObj"))

setMethod("evalGelObj", signature("momentModel", "numeric", "numeric"),
          function(object, theta, lambda, gelType, rhoFct=NULL, ...)
              {
                  gt <- evalMoment(object, theta)
                  k <- object@sSpec@k
                  if (is.null(rhoFct))
                      rhoFct <- get(paste("rho",gelType,sep=""))
                  rho <- rhoFct(gmat=gt, lambda=lambda, derive = 0, k = k[1]/k[2])
                  2*sum(rho)*k[2]/(k[1]^2*object@sSpec@bw)
              })


######################### solveGel  #########################

setGeneric("solveGel", function(object, ...) standardGeneric("solveGel"))

setMethod("solveGel", signature("momentModel"),
          function(object, gelType="EL", theta0=NULL, lambda0=NULL, lamSlv=NULL,
                   coefSlv=c("optim","nlminb","constrOptim"), rhoFct=NULL, 
                   lControl=list(), tControl=list())
          {
              coefSlv <- match.arg(coefSlv)
              if ("restrictedLam" %in% names(lControl))
              {
                  .restrictedLam <- lControl[["restrictedLam"]]
                  lControl["restrictedLam"] <- NULL
              } else {
                  .restrictedLam <- integer()
              }
              f <- function(theta, model, lambda0, slv, lcont, gelType, rhoFct,
                            returnL=FALSE, restrictedLam)
              {
                  names(theta) <- modelDims(model)$parNames
                  gt <- evalMoment(model, theta)
                  k <- model@sSpec@k
                  args <- c(list(gmat=gt, lambda0=lambda0, gelType=gelType,
                                 rhoFct=rhoFct, restrictedLam=restrictedLam),
                            lcont, k=k[1]/k[2])
                  res <- do.call(slv, args)
                  if (returnL)
                      return(res)
                  res$obj
              }
              if (is.null(lambda0))
                  lambda0 <- rep(0, modelDims(object)$q)
              if (is.null(theta0))
              {
                  if (!("theta0"%in%slotNames(object)))
                      stop("theta0 must be provided")
                  theta0 <- modelDims(object)$theta0
              }
              if (is.null(lamSlv))
                  lamSlv <- getLambda
              if (coefSlv == "nlminb")
              {
                  args <- c(list(start=theta0, objective=f, gelType=gelType,
                                 model=object, lambda0=lambda0, rhoFct=rhoFct,
                                 slv=lamSlv, lcont=lControl,
                                 restrictedLam=.restrictedLam), tControl)
              } else {
                  args <- c(list(par=theta0, fn=f, model=object, lambda0=lambda0,
                                 slv=lamSlv, lcont=lControl, gelType=gelType,
                                 rhoFct=rhoFct, restrictedLam=.restrictedLam), tControl)
              }
              res <- do.call(get(coefSlv), args)
              resl <- f(res$par,  object, lambda0, lamSlv, gelType=gelType,
                        rhoFct=rhoFct, lControl, TRUE, .restrictedLam)
              names(resl$lambda) <- modelDims(object)$momNames
              theta <- res$par
              names(theta) <- modelDims(object)$parNames                  
              list(theta=theta, convergence=res$convergence,
                   lambda=resl$lambda, lconvergence=resl$convergence,
                   restrictedLam=.restrictedLam)
          })

################ model fit ###########

setGeneric("gelFit", function(model, ...) standardGeneric("gelFit"))

setMethod("gelFit", signature("momentModel"), valueClass="gelfit", 
          definition = function(model, gelType="EL", rhoFct=NULL,
              initTheta=c("gmm", "modelTheta0"), theta0=NULL,
              lambda0=NULL, vcov=FALSE, ...)
          {
              Call <- try(match.call(call=sys.call(sys.parent())), silent=TRUE)
              if (inherits(Call,"try-error"))
                  Call <- NULL
              spec <- modelDims(model)
              initTheta = match.arg(initTheta)
              if (is.null(theta0))
              {
                  if (initTheta == "gmm")
                      theta0 <- gmmFit(model)@theta
                  else if (!is.null(spec$theta0))
                      theta0 <- spec$theta0
                  else
                      stop("starting values is missing for the coefficient vector")
              }
              res <- solveGel(model, theta0=theta0, lambda0=lambda0,
                              gelType=gelType, rhoFct=rhoFct, ...)
              if (!is.null(rhoFct))
                  gelType <- "Other"
              gelfit <- new("gelfit", theta=res$theta, convergence=res$convergence,
                            lconvergence=res$lconvergence$convergence,
                            lambda=res$lambda, call=Call,
                            gelType=list(name=gelType, rhoFct=rhoFct),
                            vcov=list(), model=model,
                            restrictedLam = res$restrictedLam)
              if (vcov)
                  gelfit@vcov <- vcov(gelfit)
              gelfit
          })


############# eval ###########

setGeneric("evalGel", function(model, ...) standardGeneric("evalGel"))

setMethod("evalGel", signature("momentModel"),
          function(model, theta, lambda=NULL, gelType="EL", rhoFct=NULL,
                   lamSlv=NULL, lControl=list(), ...)
          {
              Call <- try(match.call(call=sys.call(sys.parent())), silent=TRUE)
              if (inherits(Call,"try-error"))
                  Call <- NULL
              spec <- modelDims(model)
              theta <- setCoef(model, theta)
              type <- paste("Eval-", gelType, sep="")
              if (is.null(lambda))
              {
                  if ("restrictedLam" %in% names(lControl))
                  {
                      .restrictedLam <- lControl[["restrictedLam"]]
                      lControl["restrictedLam"] <- NULL
                  } else {
                      .restrictedLam <- integer()
                  }
                  gt <- evalMoment(model, theta)
                  k <- model@sSpec@k
                  args <- c(list(gmat=gt, gelType=gelType,
                                 rhoFct=rhoFct, restrictedLam=.restrictedLam),
                            lControl, k=k[1]/k[2])
                  if (is.null(lamSlv))
                      lamSlv <- getLambda
                  res <- do.call(lamSlv, args)
                  lambda <- res$lambda
                  lconvergence <- res$convergence$convergence
                  type <- paste(type, " with optimal lambda", sep="")
              } else {
                  lconvergence <- 1
                  type <- paste(type, " with fixed lambda", sep="")
                  .restrictedLam <- integer()
              }
              names(lambda) <- spec$momNames
              if (!is.null(rhoFct))
                  gelType <- "Other"
              new("gelfit", theta=theta, convergence=1, lconvergence=lconvergence,
                  lambda=lambda, call=Call, gelType=list(name=gelType, rhoFct=rhoFct),
                  vcov=list(), model=model, restrictedLam = .restrictedLam)
          })


