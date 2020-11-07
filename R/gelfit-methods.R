####  All methods for gelfit class
#####################################

### Hidden functions

.minvTest <- function (object, which = 1:2, fact = 2, npoints=30, level=0.95,
                       type=c("LR", "LM", "J"), cores=8) 
{
    type <- match.arg(type)
    if (length(which) != 2) 
        stop("You must select 2 coefficients")
    if (length(coef(object)) < 2) 
        stop("You need at least two coefficients")
    W <- confint(object, parm=which, level=level, area=TRUE, npoints=npoints,
                 fact=fact, type="Wald")
    p <- W@areaPoints
    spec <- modelDims(object@model)
    df <- spec$q-spec$k
    if (df > 0)
    {
        test0 <- c(specTest(object, type=type)@test)
        test0 <- test0[1]
    } else {
        test0 <- 0
    }
    f <- function(delta, pti, obj, which, type, test0, level)
    {
        b <- coef(obj)[which]
        pti <- b*(1-delta) + pti*delta
        R <- paste(names(b), "=", pti, sep="")
        if (obj@call[[1]] != "gelFit")
        {
            fit <- suppressWarnings(update(obj, cstLHS=R))
        } else {
            mod <- restModel(obj@model, R)
            fit <- suppressWarnings(update(obj, newModel=mod))
        }
        test <- c(specTest(fit, type=type)@test)[1]-test0
        test-qchisq(level, 2)
    }
    res <- try(mclapply(1:nrow(p), function(i) {
        r <- try(uniroot(f, c(0,fact), pti=p[i,], obj=object, which=which, type=type,
                         test0=test0, level=level), silent=TRUE)
        b <- coef(object)[which]
        if (inherits(r, "try-error"))
            c(NA,NA)
        else
            b*(1-r$root) + p[i,]*r$root        
    }, mc.cores=cores))
    do.call(rbind, res)
}

.invTest <- function(object, which, level = 0.95, fact = 3,
                     type=c("LR", "LM", "J"), corr=NULL, ...)
{
    type <- match.arg(type)
    if (length(which) > 1)
        stop("tests are inverted only for one parameter")
    spec <- modelDims(object@model)
    df <- spec$q-spec$k
    if (df > 0)
    {
        test0 <- c(specTest(object, type=type, ...)@test)
        test0 <- test0[1]
    } else {
        test0 <- 0
    }
    v <- diag(vcov(object, ...)$vcov_par)
    sdcoef <- sqrt(v[which])
    coef <- coef(object)[which]
    int1 <- c(coef, coef + fact*sdcoef)
    int2 <- c(coef - fact*sdcoef, coef)
    fct <- function(coef, which, type, fit, level, test0, corr=NULL)
    {
        spec <- modelDims(fit@model)
        ncoef <- spec$parNames[which]
        R <- paste(ncoef, "=", coef)
        if (fit@call[[1]] != "gelFit")           
        {
            fit2 <- suppressWarnings(update(fit, cstLHS=R))
        } else {
            model <- restModel(fit@model, R)
            fit2 <- suppressWarnings(update(fit, newModel=model,
                                            theta0=coef(fit)[-which]))
        }
        test <- specTest(fit2, type=type, ...)@test[1] - test0
         if (is.null(corr))
            level - pchisq(test, 1)
        else
            level - pchisq(test/corr, 1)
    }
    res1 <- try(uniroot(fct, int1, which = which, type=type, level=level,
                        fit=object, test0=test0, corr=corr),
                silent=TRUE)
    res2 <- try(uniroot(fct, int2, which = which, type=type, level=level,
                        fit=object, test0=test0, corr=corr),
                silent=TRUE)
    if (any(c(class(res1), class(res2)) == "try-error"))
    {
        test <- c(NA,NA)
        mess <- "Could not compute the confidence interval because: \n"
        if (inherits(res1,"try-error"))
            mess <- paste(mess, "(1) ", res1[1], "\n", sep="")
        if (inherits(res2,"try-error"))
            mess <- paste(mess, "(2) ", res2[1], "\n", sep="")
        warning(mess)        
    } else {
        test <- sort(c(res1$root, res2$root))
    }
    test
}
        
## coef

setMethod("coef", "gelfit", function(object) object@theta)

## print

setMethod("print", "gelfit",
          function(x, model=TRUE, lambda=TRUE, ...) {
              theta <- coef(x)
              if (model)
                  print(x@model)
              type <- x@gelType$name
              spec <- modelDims(x@model)
              if (spec$q==spec$k && type != "eval")
                  type <- paste("Just-Identified ", type, sep="")
              cat("\nEstimation: ", type,"\n")
              cat("Convergence Theta: ", x@convergence, "\n")
              cat("Convergence Lambda: ", x@lconvergence, "\n")
              if (length(x@restrictedLam))
              {
                  cat("Lambda's fixed at 0: ",
                      paste(x@restrictedLam, collapse=", ", sep=""),
                      "\n", sep="")
              }
              cat("coefficients:\n")
              print.default(format(theta, ...), print.gap=2L, quote=FALSE)
              if (lambda)
              {
                  cat("lambdas:\n")
                  print.default(format(x@lambda, ...), print.gap=2L, quote=FALSE)
              }
          })

## show

setMethod("show","gelfit", function(object) print(object))

## residuals

setMethod("residuals", "gelfit", function(object) {
    residuals(object@model, object@theta)})

## getImpProb

setGeneric("getImpProb", function(object, ...) standardGeneric("getImpProb"))

setMethod("getImpProb", "gelfit",
          function(object, posProb=FALSE, normalize=TRUE) {
              type <- object@gelType              
              if (is.null(type$rhoFct))
              {
                  if (type$name %in% c("ETEL","ETHD"))
                      type$name <- "ET"
                  rhoFct <- get(paste("rho", type$name, sep=""))
              } else {
                  rhoFct <- type$rhoFct
              }
              gt <- evalMoment(object@model, object@theta)
              k <- object@model@sSpec@k
              pt <- -rhoFct(gt, object@lambda, 1, k[1]/k[2])/nrow(gt)
              if (type$name == "EEL"  && posProb) {
                  eps <- -length(pt) * min(min(pt), 0)
                  pt <- (pt + eps/length(pt))/(1 + eps)
              }
              if (normalize)
                  pt <- pt/sum(pt)
              convMom <- colSums(pt * gt)
              convProb <- abs(sum(as.numeric(pt))-1)
              list(pt=pt, convMom=convMom, convProb=convProb)
          })

### To be removed once the above has need tested enough
#setMethod("getImpProb", "gelfit",
#          function(object) {
#              rhoFct <- object@model@gelType
#              if (is.null(rhoFct$fct))
#                  rhoFct <- get(paste("rho", rhoFct$name, sep=""))
#              else
#                  rhoFct <- rhoFct$fct
#              gt <- evalMoment(object@model, object@theta)
#              k <- object@model@wSpec$k
#              pt <- -rhoFct(gt, object@lambda, 1, k[1]/k[2])/nrow(gt)
#              if (object@model@gelType$name == "EEL") {
#                  eps <- -length(pt) * min(min(pt), 0)
#                  pt <- (pt + eps/length(pt))/(1 + eps)
#              }
#              convMom <- colSums(pt * gt)
#              convProb <- abs(sum(as.numeric(pt))-1)
#              pt <- pt/sum(pt)
#              list(pt=pt, convMom=convMom, convProb=convProb)
#          })

## vcov

setMethod("vcov", "gelfit",
          function(object, withImpProb=FALSE, tol=1e-10, robToMiss=FALSE) {
              spec <- modelDims(object@model)
              if (robToMiss)
              {
                  eta <- c(coef(object), object@lambda)
                  names(eta) <- NULL
                  mod <- momentModel(g=momFct, x=object, theta0=eta, vcov="MDS")
                  fit <- evalGmm(mod, theta=eta)
                  v <- vcov(fit)
                  spec <- modelDims(object@model)
                  Sigma <- v[1:spec$k, 1:spec$k]
                  dimnames(Sigma) <- list(spec$parNames, spec$parNames)
                  SigmaLam <- v[(spec$k+1):nrow(v), (spec$k+1):ncol(v)]
                  dimnames(SigmaLam) <- list(spec$momNames, spec$momNames)
                  return(list(vcov_par = Sigma, vcov_lambda = SigmaLam))
              }              
              q <- spec$q
              gt <- evalMoment(object@model, object@theta)
              n <- nrow(gt)
              bw <- object@model@sSpec@bw
              k <- object@model@sSpec@k
              if (withImpProb)
              {
                  pt <- getImpProb(object)$pt
                  G <- evalDMoment(object@model, object@theta, pt)
                  G <- G/k[1]
                  gt <- gt * sqrt(pt * bw/k[2])
              } else {
                  G <- evalDMoment(object@model, object@theta)
                  G <- G/k[1]
                  gt <- gt * sqrt(bw/k[2]/n)
              }
              qrGt <- qr(gt)
              piv <- qrGt$pivot
              R <- qr.R(qrGt)
              X <- forwardsolve(t(R), G[piv,])
              Y <- forwardsolve(t(R), diag(q)[piv,])
              res <- lm.fit(as.matrix(X), Y)
              u <- res$residuals
              Sigma <- chol2inv(res$qr$qr)/n
              diag(Sigma)[diag(Sigma) < 0] <- tol
              if (q == ncol(G)) {
                  SigmaLam <- matrix(0, q, q)
              } else {
                  SigmaLam <- crossprod(Y, u)/n * bw^2
                  diag(SigmaLam)[diag(SigmaLam) < 0] <- tol
              }
              piv <- sort.int(piv, index.return = TRUE)$ix
              list(vcov_par = Sigma, vcov_lambda = SigmaLam)
          })

## Summary


setMethod("summary","gelfit",
          function (object, ...) 
              {
                  if (length(object@vcov) == 0)
                      v <- vcov(object, ...)
                  else
                      v <- object@vcov
                  se.t <- sqrt(diag(v$vcov_par))
                  se.l <- sqrt(diag(v$vcov_lambda))
                  theta <- object@theta
                  lambda <- object@lambda
                  tval.t <- theta/se.t
                  tval.l <- lambda/se.l
                  coef <- cbind(theta, se.t, tval.t,
                                2*pnorm(abs(tval.t), lower.tail = FALSE))                 
                  coefl <- cbind(lambda, se.l, tval.l,
                                 2*pnorm(abs(tval.l), lower.tail = FALSE))
                  if (length(object@restrictedLam))
                      coefl[object@restrictedLam,-1] <- NA
                  stest <- specTest(object)
                  dimnames(coef) <- list(names(theta), c("Estimate", "Std. Error", 
                                                         "t value", "Pr(>|t|)"))
                  dimnames(coefl) <- list(names(lambda), c("Estimate", "Std. Error", 
                                                           "t value", "Pr(>|t|)"))
                  pt <- getImpProb(object)
                      
                  ans <- new("summaryGel", coef = coef, specTest = stest,
                             model = object@model, lambda=coefl,
                             convergence=object@convergence,gelType=object@gelType,
                             lconvergence=object@lconvergence, impProb=pt,
                             restrictedLam=object@restrictedLam)
                  ans})

## confint

setMethod("confint", "gelfit",
          function (object, parm, level = 0.95, lambda = FALSE,
                    type = c("Wald", "invLR", "invLM", "invJ"),
                    fact = 3, corr = NULL, vcov=NULL,
                    area = FALSE, npoints = 20, cores = 4, ...) 
          {
              type <- match.arg(type)
              if(.Platform$OS.type == "windows")
                  cores <- 1L
              spec <- modelDims(object@model)
              n <- spec$n
              theta <- coef(object)
              if (lambda)
              {
                  lam <- object@lambda
                  if (missing(parm))
                      parm <- 1:length(lam)                      
                  nlam <- spec$momNames
                  if (is.character(parm))
                      parm <- sort(unique(match(parm, nlam)))
                  nlam <- nlam[parm]
                  if (length(theta) == length(lam))
                  {
                      ntest <- paste("No confidence intervals for lambda",
                                     "when the model is just identified.")
                      ans <- matrix(NA, length(nlam), 2)
                  } else {
                      ntest <- "Wald confidence interval for Lambda"
                      if (is.null(vcov))
                          v <- vcov(object, ...)$vcov_lambda
                      se <- sqrt(diag(v))
                      if (missing(parm))
                          parm <- 1:length(lam)
                      se <- se[parm]
                      lam <- lam[parm]                  
                      zs <- qnorm((1 - level)/2, lower.tail = FALSE)              
                      ch <- zs * se
                      ans <- cbind(lam-ch, lam+ch)
                  }
                  dimnames(ans) <- list(nlam,
                                        c((1 - level)/2, 0.5 + level/2))
                  return(new("confint", interval=ans,
                             type=ntest, level=level, theta=lam[parm]))
              }
              if (type == "Wald")
              {
                  if (is.null(vcov))
                      v <-  vcov(object, ...)
                  return(getMethod("confint", "gmmfit")(object, parm, level,
                      vcov=v$vcov_par, area=area, npoints=npoints))
              } else {                  
                  if (missing(parm)) 
                      parm <- 1:length(theta)
                  ntheta <- spec$parNames
                  if (is.character(parm))
                      parm <- sort(unique(match(parm, ntheta)))
                  ntheta <- ntheta[parm]
                  type <- strsplit(type, "v")[[1]][2]
                  if (!area)
                  {
                      ntest <- paste("Confidence interval based on the inversion of the ", 
                                     type, " test", sep = "")
                      ans <- lapply(parm, function(w)
                          .invTest(object, w, level = level, 
                                   fact = fact, type = type, corr = corr, ...))
                      ans <- do.call(rbind, ans)
                      dimnames(ans) <- list(ntheta, c((1 - level)/2, 0.5 + level/2))
                  } else {
                      ntest <-  paste(type, " type confidence region", sep="")
                      ans <- .minvTest(object, parm, fact, npoints, level, type, cores)
                  }
              }
              if (!area)
                  new("confint", interval=ans, type=ntest, level=level, theta=theta[parm])
              else
                  new("mconfint", areaPoints=ans, type=ntest, level=level, theta=theta[parm])
          })
 
setMethod("confint", "numeric",
          function (object, parm, level = 0.95, gelType="EL", 
                    type = c("Wald", "invLR", "invLM", "invJ"),
                    fact = 3, vcov="iid")
          {
              Call <- try(match.call(call=sys.call(sys.parent())), silent=TRUE)
              if (inherits(Call,"try-error"))
                  Call <- NULL
              type <- match.arg(type)
              object <- as.data.frame(object)
              if (!is.null(Call))
                  names(object) <- as.character(Call)[2]
              g <- as.formula(paste(names(object),"~1",sep=""))
              n <- nrow(object)
              s <- sd(object[[1]], na.rm=TRUE)/sqrt(n)
              m <- mean(object[[1]], na.rm=TRUE)
              mod <- momentModel(g, ~1, vcov=vcov, data=object)
              fit <- gelFit(model=mod, gelType=gelType,
                            tControl=list(method="Brent",lower=m-s,upper=m+s))
              ans <- confint(fit, parm=1, level=level, type=type, fact=fact)
              rownames(ans@interval) <- names(object)
              ans
          })

setMethod("confint", "data.frame",
          function (object, parm, level = 0.95, gelType="EL", 
                    type = c("Wald", "invLR", "invLM", "invJ"),
                    fact = 3, vcov="iid", 
                    npoints=10, cores=4) 
          {
              type <- match.arg(type)
              if(.Platform$OS.type == "windows")
                  cores <- 1L
              if (missing(parm))
                  parm <- 1
              if (length(parm) == 1)
              {
                  x <- object[[parm]]
                  ans <- confint(x, level=level, gelType=gelType,
                                 type=type, fact=fact, vcov=vcov)
                  rownames(ans@interval) <- names(object[parm])
                  return(ans)
              }
              if (length(parm) != 2)
                  stop("You can only select 2 variable from your data.frame")
              object <- object[parm]
              parNames <- paste("mu_", names(object), sep="")
              g <- list(as.formula(paste(names(object)[1], "~", parNames[1], sep="")),
                        as.formula(paste(names(object)[2], "~", parNames[2], sep="")))
              theta0 <- colMeans(object)
              names(theta0) <- parNames
              mod <- momentModel(g, vcov=vcov, data=object,
                                 theta0=theta0)
              fit <- gelFit(mod, gelType=gelType)
              confint(fit, parm=1:2, level=level, lambda=FALSE,
                      type=type, fact=fact, corr=NULL, vcov=NULL, area=TRUE,
                      npoints=npoints, cores=cores)
          })

setMethod("confint", "matrix",
          function(object, parm, level = 0.95, gelType="EL", 
                    type = c("Wald", "invLR", "invLM", "invJ"),
                    fact = 3, vcov="iid", 
                    npoints=10, cores=4)
          {
              object <- as.data.frame(object)
              type <- match.arg(type)
              confint(object, parm, level, gelType, type, fact, vcov,
                      npoints, cores)
          })
                   

setMethod("confint", "ANY",
          function (object, parm, level = 0.95, ...) 
          {
              stats::confint(object, parm, level, ...)
          })


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


## specTest

           
setMethod("specTest", signature("gelfit", "missing"),
          function(object, which, type=c("All", "LR", "LM", "J"))
          {
              type <- match.arg(type)
              spec <- modelDims(object@model)
              gelType <- object@gelType
              q <- spec$q-length(object@restrictedLam)
              n <- spec$n
              df <- q-spec$k              
              test <- numeric()
              if (type %in% c("All","LR"))
              {
                  LR <- evalGelObj(object@model, theta=object@theta, lambda=object@lambda,
                                   gelType=gelType$name, rhoFct=gelType$rhoFct)
                  test <- c(test, LR)
                  names(test) <- "LR: "
              }
              if (type %in% c("All","LM","J"))
                  gt <- evalMoment(object@model, object@theta)
              if (type %in% c("All","LM"))
              {
                  kHat <- crossprod(gt)/n
                  LM <- n * crossprod(object@lambda, crossprod(kHat, object@lambda))/
                      (object@model@sSpec@bw^2)
                  test <- c(test, LM)
                  names(test)[length(test)] <- "LM: "                  
              }
              if (type %in% c("All","J"))
              {
                  J <- sum(lm.fit(gt, rep(1,n))$fitted.values)
                  test <- c(test, J)
                  names(test)[length(test)] <- " J: "                  
              }
              if (df == 0)
                  pv <- NA
              else
                  pv <- 1-pchisq(test, df)
              test <- cbind(test, df, pv)
              colnames(test) <- c("Statistics", "df", "pvalue")
              new("specTest", test=test, testname="Test E(g)=0")
          })

setMethod("print", "summaryGel",
          function(x, digits=5, lambda=TRUE, ...)
          {
              print(x@model)
              type <- x@gelType$name
              spec <- modelDims(x@model)              
              if (spec$q==spec$k && type != "eval")
                  type <- paste("Just-Identified ", type, sep="")
              cat("\nEstimation: ", type,"\n")            
              cat("Convergence Theta: ", x@convergence, "\n", sep="")
              cat("Convergence Lambda: ", x@lconvergence, "\n", sep="")
              if (length(x@restrictedLam))
              {
                  cat("Lambda's fixed at 0: ",
                      paste(x@restrictedLam, collapse=", ", sep=""),
                      "\n", sep="")
              }              
              cat("Average |Sum of pt*gt()]|: ", format(mean(abs(x@impProb$convMom)),
                                                        digits=5), "\n", sep="")
              cat("|Sum of pt - 1|: ", format(mean(abs(x@impProb$convProb)),
                                              digits=5), "\n", sep="")
              
              cat("\ncoefficients:\n")
              printCoefmat(x@coef, digits=digits, ...)
              if (lambda)
                  {
                      cat("\nLambdas:\n")
                      printCoefmat(x@lambda, digits=digits, ...)
                  }
              print(x@specTest)
          })


## show
setMethod("show", "summaryGel", function(object) print(object)) 
    

## update    

setMethod("update", "gelfit",
          function(object, newModel=NULL, ..., evaluate=TRUE)
          {
              if (is.null(call <- getCall(object)))
                  stop("No call argument")
              if (call[[1]] != "gelFit")
                  return(stats::update(object, ..., evaluate=evaluate))
              if (!is.null(newModel))
                  return(stats::update(object, model=newModel, ..., evaluate=evaluate))
              arg <- list(...)
              ev <- new.env(parent.frame())
              theta0 <- arg$theta0
              model <- if(is.null(newModel))
                           object@model
                       else
                           newModel
              model <- update(model, ...)
              ev[["model"]] <- model
              call[["model"]] <- quote(model)
              arg <- arg[which(is.na(match(names(arg),
                                           c("rhoFct", slotNames(model)))))]
              spec <- modelDims(model)
              if (!is.null(call[["theta0"]]))
              {
                  call[["theta0"]] <- if (is.null(theta0))
                                          spec$theta0
                                      else
                                          theta0
              } else if (!is.null(theta0)) {
                  call[["theta0"]] <- theta0
              }
              if (length(arg) > 0) 
                  for (n in names(arg)) call[[n]] <- arg[[n]]
              if (evaluate)
                  eval(call, ev)
              else
                  call
          })


### This method is for specific moment functions

setGeneric("momFct", function(eta, object, ...) standardGeneric("momFct"))

## That moment function is used to rewrite GEL models into
## GMM just identified models. It is useful to compute robust-to-misspecification s.e.

setMethod("momFct", signature("numeric", "gelfit"),
          function(eta, object) {
              spec <- modelDims(object@model)
              if (length(eta) != (spec$k+spec$q))
                  stop("eta must include theta and lambda")
              object@theta <- head(eta, spec$k)
              names(object@theta) <- spec$parNames
              object@lambda <- tail(eta, spec$q)
              names(object@lambda) <- spec$momNames
              pt <- getImpProb(object, FALSE, FALSE)$pt*spec$n
              gt <- evalMoment(object@model, object@theta)*pt
              Gtl <- evalDMoment(object@model, object@theta, pt, object@lambda)
              cbind(Gtl, gt)
          })


