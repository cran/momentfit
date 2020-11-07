rhoEL <- function(gmat, lambda, derive = 0, k = 1) 
    {
        lambda <- c(lambda)*k
        gmat <- as.matrix(gmat)
        gml <- c(gmat %*% lambda)
        switch(derive+1,
               log(1 - gml),
               -1/(1 - gml),
               -1/(1 - gml)^2)               
    }

rhoET <- function(gmat, lambda, derive = 0, k = 1) 
    {
        lambda <- c(lambda)*k
        gmat <- as.matrix(gmat)
        gml <- c(gmat %*% lambda)
        switch(derive+1,
               -exp(gml)+1,
               -exp(gml),
               -exp(gml))               
    }


rhoETEL <- function(gmat, lambda, derive = 0, k = 1) 
{
    lambda <- c(lambda)*k
    gmat <- as.matrix(gmat)
    gml <- c(gmat %*% lambda)
    w <- -exp(gml)
    w <- w/sum(w)
    n <- nrow(gmat)
    switch(derive+1,
           -log(w*n),
           NULL,
           NULL)               
}

rhoETHD <- function(gmat, lambda, derive = 0, k = 1) 
{
    lambda <- c(lambda)*k
    gmat <- as.matrix(gmat)
    gml <- c(gmat %*% lambda)
    w <- -exp(gml)
    w <- w/sum(w)
    n <- nrow(gmat)
    switch(derive+1,
    (sqrt(w)-1/sqrt(n))^2,
    NULL,
    NULL)               
}

rhoEEL <- function(gmat, lambda, derive = 0, k = 1) 
    {
        lambda <- c(lambda)*k
        gmat <- as.matrix(gmat)
        gml <- c(gmat %*% lambda)
        switch(derive+1,
               -gml - 0.5 * gml^2,
               -1 - gml,
               rep(-1, nrow(gmat)))               
    }

rhoREEL <- function(gmat, lambda, derive = 0, k = 1)
{
    rhoEEL(gmat, lambda, derive, k)
}

rhoHD <- function(gmat, lambda, derive = 0, k = 1) 
    {
        lambda <- c(lambda)*k
        gmat <- as.matrix(gmat)
        gml <- c(gmat %*% lambda)
        switch(derive+1,
               -1/(1 + gml)+1,
               1/((1 + gml)^2),
               -2/((1 + gml)^3))               
    }

Wu_lam <- function(gmat, tol=1e-8, maxiter=50, k=1)
    {
        gmat <- as.matrix(gmat)
        res <- .Fortran(F_wu, as.double(gmat), as.double(tol),
                        as.integer(maxiter), as.integer(nrow(gmat)),
                        as.integer(ncol(gmat)), as.double(k),
                        conv=integer(1), obj=double(1),
                        lambda=double(ncol(gmat)))
        list(lambda=res$lambda, convergence=list(convergence=res$conv),
             obj = res$obj)
    }

REEL_lam <- function(gmat, tol=NULL, maxiter=50, k=1)
{
        gmat <- as.matrix(gmat)
        n <- nrow(gmat)
        q <- ncol(gmat)
        res <- try(.Fortran(F_lamcuep, as.double(gmat),
                            as.integer(n), as.integer(q), as.double(k),
                            as.integer(maxiter),conv=integer(1),
                            lam=double(q),pt=double(n),
                            obj=double(1)
                            ), silent=TRUE)
        if (inherits(res,"try-error"))
            return(list(lambda=rep(0,q), obj=0, pt=rep(1/n,n),
                        convergence=list(convergence=3)))
        list(lambda=res$lam, obj=res$obj, pt=res$pt,
             convergence=list(convergence=res$conv))
    }

EEL_lam <- function(gmat, k=1)
{
    q <- qr(gmat)
    n <- nrow(gmat)
    lambda0 <- -qr.coef(q, rep(1,n))
    conv <- list(convergence=0)
    list(lambda = lambda0, convergence = conv,
         obj =  mean(rhoEEL(gmat,lambda0,0,k)))
}

ETXX_lam <- function(gmat, lambda0, k, gelType, algo, method, control)
{
    res <- getLambda(gmat, lambda0=lambda0, gelType="ET", algo=algo,
                     control=control, method=method, k=k)
    rhoFct <- get(paste("rho",gelType,sep=""))
    res$obj <- mean(rhoFct(gmat, res$lambda, 0, k))
    res
}

getLambda <- function (gmat, lambda0=NULL, gelType=NULL, rhoFct=NULL, 
                       tol = 1e-07, maxiter = 100, k = 1, method="BFGS", 
                       algo = c("nlminb", "optim", "Wu"), control = list(),
                       restrictedLam=integer()) 
{
    if (!is.null(gelType))
    {
        if (length(algo) == 3 & gelType == "EL")
            algo <- "Wu"
    }
    algo <- match.arg(algo)
    gmat <- as.matrix(gmat)    
    if (length(restrictedLam))
    {
        if (length(restrictedLam) > ncol(gmat))
            stop("The number of restricted Lambda exceeds the number of moments")
        if (!all(restrictedLam %in% (1:ncol(gmat))))
            stop(paste("restrictedLam must be a vector of integers between 1 and ",
                       ncol(gmat), sep=""))
        gmat <- gmat[,-restrictedLam,drop=FALSE]
    } else {
        restrictedLam <- integer()
    }
    if (is.null(lambda0))
    {
        lambda0 <- rep(0, ncol(gmat))
    } else {
        if (length(restrictedLam))
            lambda0 <- lambda0[-restrictedLam]
    }
    if (is.null(rhoFct))
    {
        if (is.null(gelType))
            stop("Without rhoFct, gelType must be given")
        rhoFct <- get(paste("rho",gelType,sep=""))
    } else {
        gelType <- "Other"
    }
    if (algo == "Wu" & gelType != "EL") 
        stop("Wu (2005) algo to compute Lambda is for EL only")
    if ((algo == "Wu") | (gelType %in% c("EEL","REEL","ETEL","ETHD")))
    {
        if (algo == "Wu")
            res <- Wu_lam(gmat, tol, maxiter, k)
        if (gelType == "EEL")
            res <- EEL_lam(gmat, k)
        if (gelType == "REEL")
            res <- REEL_lam(gmat, NULL, maxiter, k)
        if (gelType %in% c("ETEL", "ETHD"))
            res <- ETXX_lam(gmat, lambda0, k, gelType, algo, method, control)
    } else {    
        fct <- function(l, X, rhoFct, k) {
            r0 <- rhoFct(X, l, derive = 0, k = k)
            -mean(r0)
        }
        Dfct <-  function(l, X, rhoFct, k)
        {
            r1 <- rhoFct(X, l, derive = 1, k = k)
            -colMeans(r1 * X)
        }
        DDfct <-  function(l, X, rhoFct, k)
        {
            r2 <- rhoFct(X, l, derive = 2, k = k)
            -crossprod(X * r2, X)/nrow(X)
        }
        if (algo == "optim") {
            if (gelType == "EL")
            {
                ci <- -rep(1, nrow(gmat))
                res <- constrOptim(lambda0, fct, Dfct, -gmat, ci, control = control,
                                   X = gmat, rhoFct = rhoFct, k = k)
            } else if (gelType == "HD") {
                ci <- -rep(1, nrow(gmat))
                res <- constrOptim(lambda0, fct, Dfct, -gmat, ci, control = control,
                                   X = gmat, rhoFct = rhoFct, k = k)
            } else {
                res <- optim(lambda0, fct, gr = Dfct, X = gmat, rhoFct = rhoFct,
                             k = k, method = method, control = control)
            }
        } else {
            res <- nlminb(lambda0, fct, gradient = Dfct, hessian = DDfct,
                          X = gmat, rhoFct = rhoFct, k = k, control = control)
        }
        lambda0 <- res$par
        if (algo == "optim") 
            conv <- list(convergence = res$convergence, counts = res$counts, 
                         message = res$message)
        else
            conv <- list(convergence = res$convergence, counts = res$evaluations, 
                         message = res$message)    
        res <- list(lambda = lambda0, convergence = conv,
                    obj= mean(rhoFct(gmat,lambda0,0,k)))
    }
    if (length(restrictedLam))
    {
        lambda <- numeric(ncol(gmat)+length(restrictedLam))
        lambda[-restrictedLam] <- res$lambda
        res$lambda <- lambda
    }
    res
}

