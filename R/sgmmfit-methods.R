####  All methods for sgmmfit class
#####################################                        

## coef

setMethod("coef", "sgmmfit",
          function(object, asvector=FALSE)
              {
                  if (asvector)
                      {
                          spec <- modelDims(object@model)
                          .tetReshape(object@theta,
                                      spec$eqnNames, spec$parNames)
                      } else {
                          object@theta
                      }
              })

## print

setMethod("print", "sgmmfit",
          function(x, ...)
          {
              print(x@model)
              theta <- coef(x)
              ntype <- matrix(c("Two-Step GMM", "Iterated GMM", "CUE", 
                                "One-Step GMM with fixed weights",
                                "Equation by Equation Two-Stage Least Squares", 
                                "Evaluated at a fixed Theta; No estimation",
                                "Equation by Equation Two-Step GMM",
                                "Equation by Equation Iterated GMM",
                                "Equation by Equation CUE",
                                "One-Step GMM with fixed weights",
                                "Three-Stage Least Squares",
                                "Full-Information Instrumental Variables Efficient",
                                "Seemingly Unrelated Regression",
                                "twostep", "iter", "cue", "onestep", "tsls", "eval",
                                "EBEtwostep", "EBEiter", "EBEcue", "EBEonestep",
                                "3SLS", "FIVE", "SUR"),
                              ncol = 2)
              type <- ntype[match(x@type, ntype[, 2]), 1]
              spec <- modelDims(x@model)
              if (all(spec$q == spec$k) && x@type != "eval") 
                  type <- "Equation by Equation One-Step: Just-Identified"
              cat("\nEstimation: ", type, "\n", sep="")
              if (!is.null(x@convergence)) 
                  cat("Convergence Optim: ", x@convergence, "\n")
              if (!is.null(x@convIter)) 
                  cat("Convergence Iteration: ", x@convIter, "\n")
              cat("coefficients:")
              for (i in 1:length(theta))
                  {
                      cat("\n", x@model@eqnNames[i], ":\n", sep="")
                      print.default(format(theta[[i]], ...), print.gap = 2L, quote = FALSE)
                  }
          })         
## show
 
setMethod("show", "sgmmfit", function(object) print(object))         

## residuals

setMethod("residuals", "sgmmfit", function(object) {
    residuals(object@model, object@theta)})

## vcov

setMethod("vcov", "sgmmfit", 
          function (object, sandwich = NULL, df.adj = FALSE, 
                    breadOnly = FALSE, modelVcov=NULL) 
              {
                  if (!is.null(modelVcov))
                      {
                          if (modelVcov != object@model@vcov)
                              {
                                  slot(object@model, "vcov") <- modelVcov
                                  sandwich <- TRUE
                              }
                      }
                  spec <- modelDims(object@model)
                  if (breadOnly) {
                      vcov <- bread(object)/object@model@n
                      attr(vcov, "type") <- list(sandwich = FALSE, df.adj = FALSE, 
                                                 breadOnly = TRUE)
                      return(vcov)
                  }
                  if (is.null(sandwich)) 
                      sandwich <- !object@efficientGmm
                  meat <- meatGmm(object, sandwich)
                  if (!sandwich) {
                      vcov <- solve(meat)/spec$n
                  } else {
                      if (all(spec$k == spec$q)) {
                          G <- evalDMoment(object@model, coef(object))
                          v <- vcov(object@model, coef(object))
                          b <- lapply(G, solve)
                          b <- .GListToMat(b)
                          vcov <- b %*% v %*% t(b)/object@model@n
                      } else {
                          b <- bread(object)
                          vcov <- b %*% meat %*% b/object@model@n
                      }
                  }
                  tn <- paste(rep(spec$eqnNames, spec$k), ".", 
                              do.call("c", spec$parNames), sep = "")
                  dimnames(vcov) <- list(tn, tn)
                  if (df.adj) 
                      vcov <- vcov * spec$n/(spec$n - sum(spec$k))
                  attr(vcov, "type") <- list(sandwich = sandwich, df.adj = df.adj, 
                                             breadOnly = breadOnly)
                  vcov
              })

## meatGmm

setMethod("meatGmm", "sgmmfit",
          function(object, robust = FALSE) 
          {
              spec <- modelDims(object@model)
              G <- evalDMoment(object@model, coef(object))
              full <- all(sapply(1:length(G), function(i) ncol(G[[i]])==sum(spec$k)))
              G <- .GListToMat(G, full)
              if (!robust) {
                  wObj <- evalWeights(object@model, coef(object), "optimal")
                  meat <- quadra(wObj, G)
              }  else {
                  wObj <- object@wObj
                  v <- vcov(object@model, coef(object))
                  T1 <- quadra(wObj, v, G)
                  meat <- quadra(wObj, G, T1)
              }
              meat
          })

## bread

setMethod("bread", "sgmmfit",
          function (x, ...) {
              G <- evalDMoment(x@model, x@theta)
              wObj <- x@wObj
              spec <- modelDims(x@model)
              if (all(spec$q == spec$k) && is.character(wObj@w))
                  {
                      G <- lapply(1:length(x@model@eqnNames), function(i)
                          crossprod(solve(G[[i]])))
                      b <- .GListToMat(G)
                  } else {
                      G <- .GListToMat(G)
                      b <- solve(quadra(wObj, G))
                  }
              b
          })

## specTest 

setMethod("specTest", c("sgmmfit","missing"),
          function(object, which, df.adj = FALSE, wObj = NULL) 
    {
        spec <- modelDims(object@model)
        J_test <- "J-Test"
        if (is.null(wObj)) 
            wObj <- evalWeights(object@model, coef(object), "optimal")
        else
            J_test <- paste(J_test, "(Weights provided by user)")
        j <- evalGmmObj(object@model, coef(object), wObj)
        if (df.adj) 
            j <- j * (spec$n - sum(spec$k))/spec$n
        df <- sum(spec$q) - sum(spec$k)
        j <- cbind(j, df)
        j <- cbind(j, ifelse(df > 0, pchisq(j, df, lower.tail = FALSE), 
                             NA))
        dimnames(j) <- list("Test E(g)=0:  ", c("Statistics", 
                                       "df", "pvalue"))
        ans <- new("specTest", test = j, testname = J_test)
        ans
    })

### summary

setMethod("summary","sgmmfit",
          function (object, testStrength=TRUE, ...) {
              spec <- modelDims(object@model)
              eqnNames <- spec$eqnNames
              neqn <- length(eqnNames)
              v <- vcov(object, ...)
              se <- sqrt(diag(v))
              par <- coef(object, TRUE)
              tval <- par/se
              pval <- 2 * pnorm(abs(tval), lower.tail = FALSE)
              se <- .tetReshape(se, eqnNames, spec$parNames)
              par <- .tetReshape(par, eqnNames, spec$parNames)
              tval <- .tetReshape(tval, eqnNames, spec$parNames)
              pval <- .tetReshape(pval, eqnNames, spec$parNames)
              coef <- lapply(1:length(se), function(i) {
                  b <- cbind(par[[i]], se[[i]], tval[[i]], pval[[i]])
                  dimnames(b) <- list(names(par[[i]]),
                                      c("Estimate", "Std. Error", "t value", "Pr(>|t|)"))
                  b})
              names(coef) <- eqnNames
              df.adj <- attr(v, "type")$df.adj
              stest <- specTest(object, df.adj = df.adj)
              vcovType <- switch(object@model@vcov, HAC = "HAC", iid = "OLS", 
                                 MDS = "HC")
              strength <- lapply(1:neqn, function(i) {
                  if (inherits(object@model, "slinearModel") & testStrength)
                      momentStrength(object@model[i], par[[i]], vcovType)
                  else
                      list(strength=NULL, mess=NULL)})             
              wSpec <- object@wObj@wSpec
              ans <- new("summarySysGmm", coef = coef, type = object@type, 
                         specTest = stest, strength = strength, model = object@model, 
                         df.adj = df.adj, niter = object@niter,
                         breadOnly = attr(v, "type")$breadOnly,
                         convergence = object@convergence, 
                         wSpec = wSpec, convIter = object@convIter,
                         sandwich = attr(v, "type")$sandwich)
              ans
          })


## update

## hypothesisTest

setMethod("hypothesisTest", signature("sgmmfit", "missing"),
          function(object.u, object.r, R, rhs=NULL, vcov=NULL, ...) {
              if (inherits(object.u@model, "rsysModel"))
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

setMethod("hypothesisTest", signature("missing", "sgmmfit"),
          function (object.u, object.r, wObj = NULL) 
              {
                  ntest <- "LM Test"
                  if (!inherits(object.r@model, "rsysModel")) 
                      stop("LR tests can only be performed on restricted models")
                  b <- coef(object.r@model, coef(object.r))
                  if (is.null(wObj)) 
                      wObj <- object.r@wObj
                  uobj <- as(object.r@model, substring(class(object.r@model), 
                                                       2))
                  G <- evalDMoment(uobj, b)
                  G <- .GListToMat(G)
                  gt <- evalMoment(object.r@model, coef(object.r))
                  gt <- do.call(cbind, gt)
                  gbar <- colMeans(gt)
                  T1 <- quadra(wObj, G, gbar)
                  T2 <- quadra(wObj, G)
                  test <- object.r@model@n * c(crossprod(T1, solve(T2, T1)))
                  df.test <- sum(modelDims(uobj)$k) - sum(modelDims(object.r@model)$k)
                  pv <- 1 - pchisq(test, df.test)
                  hypo <- getRestrict(object.r@model, b)$hypo
                  new("hypothesisTest", test = test, df = df.test, pvalue = pv, 
                      hypothesis = hypo, dist = "Chi-square", type = ntest)
              })


setMethod("hypothesisTest", signature("sgmmfit", "sgmmfit"),
          function(object.u, object.r, type=c("Wald", "LR", "LM"),
                   sameVcov=TRUE, vcov=NULL, firstStepWeight=FALSE, wObj=NULL, ...) {
              type <- match.arg(type)
              ntest <- paste(type, "Test")
              if (inherits(object.u@model, "rsysModel"))
                  stop("object.u must be a fitted unrestricted model")
              if (!inherits(object.r@model, "rsysModel"))
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


### Hausman Statistics


