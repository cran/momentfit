####  All methods for gmmfit class
#####################################
                          

## coef

setGeneric("coef")
setMethod("coef", "gmmfit", function(object) object@theta)

## print

setMethod("print", "gmmfit",
          function(x, model=TRUE, ...) {
              theta <- coef(x)
              if (model)
                  print(x@model)
              ntype <- matrix(c("Two-Step GMM", "Iterated GMM", "CUE",
                                "One-Step GMM with fixed weights","Two-Stage Least Squares",
                                "Evaluated at a fixed Theta; No estimation",
                                "One-Step Efficient M.D.E.",
                                "twostep","iter","cue","onestep","tsls", "eval","mde"),
                              ncol=2)
              type <- ntype[match(x@type, ntype[,2]),1]
              spec <- modelDims(x@model)
              if (spec$q==spec$k && x@type != "eval")
                  type <- "One-Step, Just-Identified"
              cat("\nEstimation: ", type,"\n")
              if (!is.null(x@convergence))
                  cat("Convergence Optim: ", x@convergence, "\n")              
              if (!is.null(x@convIter))
                  cat("Convergence Iteration: ", x@convIter, "\n")
              cat("coefficients:\n")
              print.default(format(theta, ...), print.gap=2L, quote=FALSE)
          })

## show

setMethod("show","gmmfit", function(object) print(object))

## residuals

setMethod("residuals", "gmmfit", function(object) {
    residuals(object@model, object@theta)})

## vcov

setMethod("vcov", "gmmfit",
          function(object, sandwich=NULL, df.adj=FALSE, breadOnly=FALSE,
                   modelVcov=NULL) {
              if (!is.null(modelVcov))
                  {
                      if (modelVcov != object@model@vcov)
                          {
                              slot(object@model, "vcov") <- modelVcov
                              sandwich <- TRUE
                          }
                  }
              spec <- modelDims(object@model)
              if (breadOnly)
                  {
                      vcov <- bread(object)/spec$n
                      attr(vcov, "type") <- list(sandwich=FALSE, df.adj=FALSE,
                                                 breadOnly=TRUE)
                      return(vcov)
                  }
              if (is.null(sandwich))
                  sandwich <- !object@efficientGmm              
              meat <- meatGmm(object, sandwich)
              if (!sandwich)
                  {
                      vcov <- solve(meat)/spec$n
                  } else {
                      if (spec$k==spec$q)
                          {
                              G <- evalDMoment(object@model, coef(object))
                              v <- vcov(object@model, coef(object))
                              b <- solve(G)
                              vcov <- b%*%v%*%t(b)/spec$n
                          } else {
                              b <- bread(object)
                              vcov <- b%*%meat%*%b/spec$n
                          }
                  }
              dimnames(vcov) <- list(spec$parNames,spec$parNames)
              if (df.adj)
                  vcov <- vcov*spec$n/(spec$n-spec$k)
              attr(vcov, "type") <- list(sandwich=sandwich, df.adj=df.adj,
                                         breadOnly=breadOnly)
              vcov
          })

## meatGmm

setGeneric("meatGmm", function(object, ...) standardGeneric("meatGmm"))

setMethod("meatGmm", "gmmfit",
          function(object, robust=FALSE) {
              G <- evalDMoment(object@model, coef(object))
              if (!robust)
              {
                  wObj <- evalWeights(object@model, coef(object), "optimal")
                  meat <- quadra(wObj, G)
              } else {
                  wObj <- object@wObj
                  v <- vcov(object@model, coef(object))                  
                  T1 <- quadra(wObj, v, G)
                  meat <- quadra(wObj, G, T1)
              }
              meat})

## bread

setGeneric("bread")
setMethod("bread", "gmmfit",  function(x, ...) {
    G <- evalDMoment(x@model, x@theta)
    wObj <-  x@wObj
    if (dim(G)[1] == dim(G)[2] && is.character(wObj@w))
        crossprod(solve(G))
    else 
        solve(quadra(wObj, G))
})


## specTest 

setGeneric("specTest", function(object, which, ...) standardGeneric("specTest"))

setMethod("specTest", signature("gmmfit", "missing"),
          function(object, which, df.adj=FALSE, wObj=NULL) {
              spec <- modelDims(object@model)
              J_test <- "J-Test"
              if (is.null(wObj))
                  wObj <- evalWeights(object@model, coef(object), "optimal")
              else
                  J_test <- paste(J_test, "(Weights provided by user)")
              j <- evalGmmObj(object@model,coef(object), wObj)
              if (df.adj)
                  j <- j*(spec$n-spec$k)/spec$n
              df <- spec$q-spec$k
              j <- cbind(j, df)
              j <- cbind(j, ifelse(df > 0, pchisq(j, df, lower.tail = FALSE), NA))
              dimnames(j) <- list("Test E(g)=0:  ", c("Statistics", "df", "pvalue"))
              ans <- new("specTest", test=j, testname=J_test)
              ans})

setMethod("specTest", signature("gmmfit", "numeric"),
          function(object, which) {
              which <- sort(unique(as.integer(which)))
              spec <- modelDims(object@model)
              q <- spec$q
              if (!all(which%in%(1:q)))
                  stop("SubMoment must be between 1 and q")
              if (spec$k > (q-length(which)))
                  stop("Th model without the tested conditions is under-identified")
              mod2 <- object@model[-which]
              w <- object@wObj
              obj2 <- gmmFit(mod2, weights=w[-which])
              J <- specTest(object, wObj=w)@test[1]
              J1 <- specTest(obj2, wObj=w[-which])@test[1]
              j <- J-J1                                        
              J_test <- paste("Testing the following subset of moments:\n\t{",
                              paste(object@model@momNames[which], collapse=","),"}",
                              sep="")
              df <- length(which)
              j <- cbind(j, df)
              j <- cbind(j, ifelse(df > 0, pchisq(j, df, lower.tail = FALSE), NA))
              dimnames(j) <- list("Test E(g)=0:  ", c("Statistics", "df", "pvalue"))
              ans <- new("specTest", test=j, testname=J_test)
              ans})

### summary

setGeneric("summary")
setMethod("summary", "gmmfit",
          function(object, testStrength=TRUE, ...)
          {
              v <- vcov(object, ...)
              se <- sqrt(diag(v))
              par <- coef(object)
              tval <- par/se
              coef <- cbind(par, se, tval, 2 * pnorm(abs(tval), 
                                                     lower.tail = FALSE))
              df.adj <- attr(v, "type")$df.adj
              stest <- specTest(object, df.adj=df.adj)
              vcovType <- switch(object@model@vcov,
                                 HAC="HAC",
                                 iid="OLS",
                                 MDS="HC",
                                 CL="CL")
              strength <-  if (testStrength){
                               momentStrength(object@model, coef(object), vcovType)
                           } else { list(strength=NULL, mess=NULL) }
              dimnames(coef) <- list(names(par), 
                                     c("Estimate", "Std. Error", "t value", "Pr(>|t|)"))
              wSpec <- object@wObj@wSpec
              ans <- new("summaryGmm", coef=coef, type=object@type,specTest=stest,
                         strength=strength, model=object@model, df.adj=df.adj,
                         niter=object@niter, breadOnly=attr(v, "type")$breadOnly,
                         convergence=object@convergence,wSpec=wSpec, 
                         convIter=object@convIter, sandwich=attr(v, "type")$sandwich)
              ans })

## confint

setGeneric("confint")

setMethod("confint","gmmfit", 
          function (object, parm, level = 0.95, vcov=NULL, area=FALSE,
                    npoints=50, ...) 
          {
              if (is.null(vcov)) 
                  vcov <- vcov(object, ...)
              theta <- coef(object)
              if (missing(parm)) 
                  parm <- 1:length(theta)
              if (is.character(parm)) 
                  parm <- sort(unique(match(parm, names(theta))))
              if (area)
              {
                  if (length(parm) != 2)
                      stop("For confidence region, length(parm) must be 2")
                  ntest <- "Wald type confidence region"
                  r <- sqrt(qchisq(level, 2))
                  v <- chol(vcov[parm,parm])
                  ang <- seq(0, 2*pi, length.out=npoints)
                  c1  <- rbind(r*cos(ang), r*sin(ang))
                  p2 <- t(crossprod(v,c1) + theta[parm])
                  colnames(p2) <- names(theta)[parm]
                  rownames(p2) <- NULL
                  obj <- new("mconfint", areaPoints=p2, type=ntest, level=level,
                             theta=theta[parm])
                  return(obj)
              }                  
              se <- sqrt(diag(vcov)[parm])
              ntest <- "Wald type confidence interval"
              theta <- theta[parm]
              zs <- qnorm((1 - level)/2, lower.tail = FALSE)
              ch <- zs * se
              ans <- cbind(theta - ch, theta + ch)
              dimnames(ans) <- list(names(theta), c((1 - level)/2, 
                                                    0.5 + level/2))
              new("confint", interval = ans, type = ntest, level = level,
                  theta=theta[parm])
          })

## confint method

setMethod("print", "confint",
          function(x, digits=4, ...)
          {
              cat("\n", x@type, "\n")
              print.default(x@interval, digits = digits, print.gap = 2, 
                            quote = FALSE)
              cat("\n")
          })

setMethod("show", "confint", function(object) print(object))


## mconfint methods

### print, show and plot for mconfint

setMethod("print", "mconfint", 
          function(x, digits=4, ...) 
          { 
              cat(x@type, "\n")
              cat(paste(rep("*", nchar(x@type)), collapse="", sep=""), "\n", sep="")
              cat("Level: ", x@level, "\n", sep="")
              cat("Number of points: ", nrow(x@areaPoints), "\n", sep="")
              cat("Variables:\n")
              v <- colnames(x@areaPoints)
              r <- formatC(apply(x@areaPoints, 2, range), digits=digits, ...)
              cat("\tRange for ", v[1], ": [", r[1,1], ", ", r[2,1], "]\n", sep="")
              cat("\tRange for ", v[2], ": [", r[1,2], ", ", r[2,2], "]\n", sep="")
              })

setMethod("show", "mconfint", function(object) print(object))

setGeneric("plot")
setMethod("plot", "mconfint", function(x, y, main=NULL, xlab=NULL, ylab=NULL, 
                                       pch=21, bg=1, Pcol=1, ylim=NULL, xlim=NULL,
                                       add=FALSE, addEstimates=TRUE, ...)
{
    v <- colnames(x@areaPoints)
    if (!add)
        {
            if (is.null(main))
                main <- paste(100*x@level, "% confidence region for ", v[1], " and ", v[2],
                              sep="")
            if (is.null(xlab))
                xlab <- v[1]
            if (is.null(ylab))
                ylab <- v[2]
            if (is.null(ylim))
                ylim <- range(x@areaPoints[,2])
            if (is.null(xlim))
                xlim <- range(x@areaPoints[,1])            
            plot(x@areaPoints, xlab=xlab, ylab=ylab, main=main, pch=pch, bg=bg,
                 ylim=ylim, xlim=xlim, col=Pcol)
            grid()
            if (addEstimates)
                {
                    points(x@theta[1], x@theta[2], pch=20)
                    text(x@theta[1], x@theta[2], expression(hat(theta)), pos=3)
                }
        } else {
            points(x@areaPoints[,1],x@areaPoints[,2],pch=pch, bg=bg, col=Pcol)
        }
    polygon(x@areaPoints[,1], x@areaPoints[,2], ...)
    invisible()
})

setMethod("plot", "ANY", function(x, y, ...) 
    if(!missing(y))
        graphics::plot(x, y, ...)
    else
        graphics::plot(x, ...)
    )


## update

setMethod("update", "gmmfit",
          function(object, ..., evaluate=TRUE)
              {
                  if (is.null(call <- getCall(object)))
                      stop("No call argument")
                  if (call[[1]] != "gmmFit")
                      return(stats::update(object, ..., evaluate=evaluate))
                  arg <- list(...)
                  if (length(arg) == 0)
                      return(object)
                  model <- update(object@model, ...)
                  ev <- new.env(parent.frame())
                  arg <- arg[which(is.na(match(names(arg), slotNames(model))))]
                  call[["model"]] <- quote(model)
                  ev$model <- model                  
                  if (length(arg) > 0)
                      for (n in names(arg))
                          call[[n]] <- arg[[n]]
                  if (evaluate)
                      eval(call, ev)
                  else
                      call
              })

## hypothesisTest

setGeneric("hypothesisTest", function(object.u, object.r, ...)
    standardGeneric("hypothesisTest"))

setMethod("hypothesisTest", signature("gmmfit", "missing"),
          function(object.u, object.r, R, rhs=NULL, vcov=NULL, ...) {
              if (inherits(object.u@model, "rmomentModel"))
                  stop("object.u must be a fitted unrestricted model")
              rest <- getRestrict(object.u@model, coef(object.u), R, rhs)
              v <- if(is.null(vcov))
                       vcov(object.u, ...)
                   else
                       vcov
              v <- rest$dR%*%v%*%t(rest$dR)
              g <- rest$R-rest$q
              test <- sum(g * c(solve(v, g)))
              df.test <- nrow(rest$dR)
              pv <- 1-pchisq(test,df.test)
              type <- "Wald Test"
              new("hypothesisTest", test=test, df=df.test, pvalue=pv, hypothesis=rest$hypo,
                  dist="Chi-square", type=type)              
          })

setMethod("hypothesisTest", signature("gmmfit", "gmmfit"),
          function(object.u, object.r, type=c("Wald", "LR", "LM"),
                   sameVcov=TRUE, vcov=NULL, firstStepWeight=FALSE, wObj=NULL, ...) {
              type <- match.arg(type)
              ntest <- paste(type, "Test")
              if (inherits(object.u@model, "rmomentModel"))
                  stop("object.u must be a fitted unrestricted model")
              if (!inherits(object.r@model, "rmomentModel"))
                  stop("object.u must be a fitted restricted model")
              rest <- getRestrict(object.r@model, coef(object.u))
              if (type == "Wald")
                  return(hypothesisTest(object.u=object.u, R=rest$orig.R, rhs=rest$orig.rhs,
                                        vcov=vcov, ...))
              if (type == "LM")
                  return(hypothesisTest(object.r=object.r, wObj=wObj))

              if (is.null(wObj))
                  {
                      wObj.r <- NULL
                      if (!all(c(object.r@efficientGmm, object.u@efficientGmm)))
                          stop("LR tests can only be performed on efficient GMM fits")
                      if (firstStepWeight)
                          {
                              wObj <- object.u@wObj
                              if (sameVcov)
                                  wObj.r <- wObj
                              else
                                  wObj.r <- object.r@wObj
                          }  else {
                              if (sameVcov)
                                  wObj.r <- evalWeights(object.u@model, coef(object.u),
                                                        "optimal")
                          }
                  } else {
                      wObj.r <- wObj
                      ntest <- paste(ntest, "(Test based on a provided weighting matrix)")
                  }
              Tu <- specTest(object.u, wObj=wObj, ...)@test[1]
              Tr <- specTest(object.r, wObj=wObj.r, ...)@test[1]
              test <- Tr-Tu
              if (test<0)
                  {
                      test <- 0
                      warning("The statistics is negative. Set sameVcov to TRUE to avoid it")
                  }
              df.test <- nrow(rest$dR)
              pv <- 1-pchisq(test,df.test)
              hypo <- rest$hypo
              new("hypothesisTest", test=test, df=df.test, pvalue=pv, hypothesis=hypo,
                  dist="Chi-square", type=ntest)              
          })

setMethod("hypothesisTest", signature("missing", "gmmfit"),
          function(object.u, object.r, wObj=NULL) {
              ntest <- "LM Test"
              if (!inherits(object.r@model, "rmomentModel"))
                  stop("LR tests can only be performed on restricted models")
              b <- coef(object.r@model, coef(object.r))
              if (is.null(wObj))
                  wObj <- object.r@wObj
              uobj <- as(object.r@model, substring(class(object.r@model), 2))
              G <- evalDMoment(uobj, b)
              gbar <- colMeans(evalMoment(object.r@model, coef(object.r)))
              T1 <- quadra(wObj, G, gbar)
              T2 <- quadra(wObj, G)
              spec <- modelDims(object.r@model)
              test <- spec$n*c(crossprod(T1, solve(T2, T1)))
              df.test <- length(b)-spec$k
              pv <- 1-pchisq(test,df.test)
              hypo <- getRestrict(object.r@model, b)$hypo
              new("hypothesisTest", test=test, df=df.test, pvalue=pv, hypothesis=hypo,
                  dist="Chi-square", type=ntest)              
          })


### Hausman Statistics


# First a solve for non invertible matrix using the Moore-Penrose generalized inverse

.solveMP <- function(A, b, tol=sqrt(.Machine$double.eps))
    {
        res <- svd(A)
        d <- res$d
        if (all(d>tol*d[1]))
            {
                r <- length(d)
                sol <- solve(A,b)
            } else {
                di <- ifelse(d<tol*d[1], 0, 1/d)
                r <- sum(d>tol*d[1])
                sol <- (res$v%*%diag(di)%*%t(res$u))%*%b
            }
        list(sol=c(sol), rank=r)
    }


setGeneric("DWH", function(object1, object2, ...) standardGeneric("DWH"))

setMethod("DWH", signature("gmmfit", "missing"),
          function(object1, object2)
              {
                  if(!inherits(object1@model, "linearModel"))
                      stop("Durbin-Wu-Hausman tests are for linear models. For Hausman test, provide two different GMM fit")
                  spec <- modelDims(object1@model)
                  EndoVars <- !(spec$parNames %in% spec$momNames)
                  vcovType <- object1@model@vcov
                  if (all(!EndoVars))
                      stop("DWH tests require at least one endogenous variable")

                  X <- model.matrix(object1@model)
                  Y <- X[,EndoVars,drop=FALSE]
                  Z <- model.matrix(object1@model, "instrument")
                  Yhat <- lm.fit(Z, Y)$fitted
                  Y <- modelResponse(object1@model)
                  res <- lm(Y~X+Yhat-1)
                  R <- diag(length(coef(res)))[-(1:ncol(X)),, drop=FALSE]
                  v <- switch(vcovType,
                              iid=vcov(res),
                              MDS=vcovHC(res,"HC1"),
                              HAC=vcovHAC(res))
                  v <- R%*%v%*%t(R)
                  T1 <- R%*%coef(res)
                  test <- t(T1)%*%solve(v, T1)
                  df.test <- nrow(R)
                  test <- cbind(test, df.test,  1-pchisq(test, df.test))
                  dimnames(test) <- list("DWH Test:  ", c("Statistics", 
                                                          "df", "pvalue"))
                  new("specTest", test=test, testname="Durbin-Wu-Hausman Test")
              })

setMethod("DWH", signature("gmmfit", "gmmfit"),
          function(object1, object2, tol=sqrt(.Machine$double.eps),
                   v1=NULL, v2=NULL, ...)
              {
                  spec1 <- modelDims(object1@model)
                  spec2 <- modelDims(object2@model)
                  if (!identical(spec1$parNames, spec2$parNames))
                      stop("You are trying to perform a DWH test on different models")
                  b1 <- coef(object1)
                  if (is.null(v1))
                      v1 <- vcov(object1, ...)
                  b2 <- coef(object2)
                  if (is.null(v2))
                      v2 <- vcov(object2, ...)
                  v <- v1-v2
                  b <- b1-b2
                  res <- .solveMP(v, b, tol)
                  test <- abs(t(b)%*%res$sol)
                  df.test <- res$rank
                  test <- cbind(test, df.test,  1-pchisq(test, df.test))
                  dimnames(test) <- list("Hausman Test:  ", c("Statistics", 
                                                          "df", "pvalue"))
                  new("specTest", test=test, testname="Hausman Test")
              })

setMethod("DWH", signature("gmmfit", "lm"),
          function(object1, object2, tol=sqrt(.Machine$double.eps),
                   v1=NULL, v2=NULL, ...)
              {
                  spec <- modelDims(object1@model)
                  if (!identical(spec$parNames, names(coef(object2))))
                      stop("You are trying to perform a DWH test on different models")
                  b1 <- coef(object1)
                  if (is.null(v1))
                      v1 <- vcov(object1, ...)
                  b2 <- coef(object2)
                  if (is.null(v2))
                      v2 <- vcov(object2)                
                  v <- v1-v2
                  b <- b1-b2
                  res <- .solveMP(v, b, tol)
                  test <- abs(t(b)%*%res$sol)
                  df.test <- res$rank
                  test <- cbind(test, df.test,  1-pchisq(test, df.test))
                  dimnames(test) <- list("Hausman Test:  ", c("Statistics", 
                                                          "df", "pvalue"))
                  new("specTest", test=test, testname="Hausman Test")
              })



