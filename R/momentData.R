######### Function to arrange the data for the gmmModel objects #################

.multiToSys <- function(formula, h, data, survOptions=list(), vcovOptions=list(),
                        na.action="na.omit")
{
    modelF <- model.frame(formula, data, na.action="na.pass",
                          drop.unused.levels=TRUE)
    Y <- model.response(modelF)
    modelF <- modelF[-1]
    Yn <- formula[[2]]
    Yn <- paste(Yn, ".", 1:ncol(Y), sep="")
    g <- lapply(1:length(Yn), function(i) {
        f <- formula
        f[[2]] <- as.symbol(Yn[i])
        f})
    colnames(Y) <- Yn
    modelF <- cbind(Y, modelF)
    if (any(class(h) == "formula"))
        {
            instF <- model.frame(h, data, na.action="na.pass",
                                 drop.unused.levels=TRUE)
        } else {
            h <- as.data.frame(h)
            chk <- apply(h, 2, function(x) all(x==x[1]))
            h <- h[, !chk]
            intercept <- any(chk)
            if (ncol(h) == 0)
                {                        
                    formula <- ~1
                } else {
                    if (is.null(colnames(h)))
                        colnames(h) <- paste("h", 1:ncol(h), sep="")
                    formh <- paste(colnames(h), collapse="+")
                    if (!intercept)
                        formh <- paste(formh, "-1", sep="")
                    formula <- as.formula(paste("~",formh))
                }
                instF <- model.frame(formula, h, na.action="na.pass",
                                     drop.unused.levels=TRUE)
        }
    h <- lapply(1:ncol(Y), function(i) formula(attr(instF, "terms"), .GlobalEnv))
    data <- cbind(modelF, instF)
    data <- data[,!duplicated(colnames(data))]
    return(.slModelData(g,h,data,survOptions, vcovOptions,na.action))
}

.lModelData <- function(formula, h, data, survOptions=list(), vcovOptions=list(),
                      na.action="na.omit")
    {
        modelF <- model.frame(formula, data, na.action="na.pass",
                              drop.unused.levels=TRUE)
        if (is.matrix(modelF[[1]]))
            return(.multiToSys(formula, h, data))
        parNames <- colnames(model.matrix(attr(modelF, "terms"), modelF))
        k <- length(parNames)
        if (any(class(h) == "formula"))
            {
                instF <- model.frame(h, data, na.action="na.pass",
                                     drop.unused.levels=TRUE)
            } else {
                h <- as.data.frame(h)
                chk <- apply(h, 2, function(x) all(x==x[1]))
                h <- h[, !chk]
                intercept <- any(chk)
                if (ncol(h) == 0)
                    {                        
                        formula <- ~1
                    } else {
                        if (is.null(colnames(h)))
                            colnames(h) <- paste("h", 1:ncol(h), sep="")
                        formh <- paste(colnames(h), collapse="+")
                        if (!intercept)
                            formh <- paste(formh, "-1", sep="")
                        formula <- as.formula(paste("~",formh))
                    }
                instF <- model.frame(formula, h, na.action="na.pass",
                                     drop.unused.levels=TRUE)
            }
        momNames <- colnames(model.matrix(attr(instF, "terms"), instF))
        q <- length(momNames)
        isEndo <- !(parNames %in% momNames)
        tmpDat <- cbind(modelF, instF)
        add  <- survOptions$weights
        if (!is.null(vcovOptions$cluster))
            add <- cbind(as.matrix(vcovOptions$cluster), add)        
        if (!is.null(add))
            tmpDat <- cbind(tmpDat, add)
        na <- attr(get(na.action)(tmpDat), "na.action")[]
        if (!is.null(na))
        {
            modelF <- modelF[-na,,drop=FALSE]
            instF <- instF[-na,,drop=FALSE]
            if (!is.null(vcovOptions$cluster))
                {
                    if (is.null(dim(vcovOptions$cluster)))
                        vcovOptions$cluster <- vcovOptions$cluster[-na]
                    else
                        vcovOptions$cluster <- vcovOptions$cluster[-na,,drop=FALSE]
                }
            if (!is.null(survOptions$weights))
                survOptions$weights <- survOptions$weights[-na]
        }
        if (is.null(na))
            na <- integer()
        n <- nrow(modelF)
        list(modelF=modelF,  instF=instF, n=n, k=k, q=q, momNames=momNames,
             parNames=parNames, isEndo=isEndo, varNames=parNames, omit=na,
             vcovOptions=vcovOptions, survOptions=survOptions)
    }

.formModelData <- function(formula, theta0, data, survOptions=list(), vcovOptions=list(),
                         na.action="na.omit")
    {
        res <- lapply(formula, function(f) .nlModelData(f, ~1, theta0, data))
        fRHS <- lapply(res, function(r) r$fRHS)
        fLHS <- lapply(res, function(r) r$fLHS)
        parNames <- res[[1]]$parNames
        varNames <- do.call("c", lapply(res, function(r) r$varNames))
        varNames <- unique(varNames)       
        chkLHS <- sapply(fLHS, function(r) any(all.vars(r) %in% names(theta0)))
        chkRHS <- sapply(fRHS, function(r) any(all.vars(r) %in% names(theta0)))
        isMDE <- all(chkLHS) |  all(chkRHS)        
        modelF <- sapply(varNames, function(n) data[[n]])
        modelF <- as.data.frame(modelF)
        tmpDat <- modelF
        add  <- survOptions$weights
        if (!is.null(vcovOptions$cluster))
            add <- cbind(as.matrix(vcovOptions$cluster), add)        
        if (!is.null(add))
            tmpDat <- cbind(tmpDat, add)
        na <- attr(get(na.action)(tmpDat), "na.action")[]
        if (!is.null(na))
        {
            modelF <- modelF[-na,,drop=FALSE]
            if (!is.null(vcovOptions$cluster))
                {
                    if (is.null(dim(vcovOptions$cluster)))
                        vcovOptions$cluster <- vcovOptions$cluster[-na]
                    else
                        vcovOptions$cluster <- vcovOptions$cluster[-na,,drop=FALSE]
                }
            if (!is.null(survOptions$weights))
                survOptions$weights <- survOptions$weights[-na]
        }
        if (is.null(na))
            na <- integer()
        k <- length(theta0)
        q <- length(formula)
        if (is.null(names(formula)))
            momNames <- paste("Mom_", 1:q, sep="")
        else
            momNames <- names(formula)
        isEndo <- rep(FALSE, length(varNames))
        n <- nrow(modelF)
        list(modelF=modelF,  fRHS=fRHS, fLHS=fLHS, n=n, k=k, q=q,
             momNames=momNames, parNames=parNames, varNames=varNames, isEndo=isEndo,
             isMDE=isMDE,omit=na, vcovOptions=vcovOptions, survOptions=survOptions)
    }



.nlModelData <- function(formula, h, theta0, data, survOptions=list(), vcovOptions=list(),
                       na.action="na.omit")
    {
        varNames <- all.vars(formula)
        parNames <- names(theta0)
        varNames <- varNames[!(varNames %in% parNames)]
        modelF <- try(sapply(varNames, function(n) data[[n]]), silent=TRUE)
        if (any(class(modelF)=="try-error"))
            stop("some variables are missing from data")
        modelF <- as.data.frame(modelF)        
        allVar <- c(as.list(modelF), as.list(theta0))
        k <- length(theta0)
        if (length(formula) == 3L)
        { 
            fLHS <- as.expression(formula[[2]])
            chk <- try(eval(fLHS, allVar))
            if (any(class(chk)=="try-error"))
                stop("Cannot evaluate the LHS")
            fRHS <- as.expression(formula[[3]])
            chk <- try(eval(fRHS, allVar))
            if (any(class(chk)=="try-error"))
                stop("Cannot evaluate the RHS")
        } else {
            fLHS <- NULL
            fRHS <- as.expression(formula[[2]])
            chk <- try(eval(fRHS, allVar))
            if (any(class(chk)=="try-error"))
                stop("Cannot evaluate the RHS")
        }
        if (any(class(h) == "formula"))
            {
                instF <- model.frame(h, data, na.action="na.pass",
                                     drop.unused.levels=TRUE)
            } else {
                h <- as.data.frame(h)
                chk <- apply(h, 2, function(x) all(x==x[1]))
                h <- h[, !chk]
                intercept <- any(chk)
                if (ncol(h) == 0)
                    {                        
                        formula <- ~1
                    } else {
                        if (is.null(colnames(h)))
                            colnames(h) <- paste("h", 1:ncol(h), sep="")
                        formh <- paste(colnames(h), collapse="+")
                        if (!intercept)
                            formh <- paste(formh, "-1", sep="")
                        formula <- as.formula(paste("~",formh))
                    }
                instF <- model.frame(formula, h, na.action="na.pass",
                                     drop.unused.levels=TRUE)
            }
        momNames <- colnames(model.matrix(attr(instF, "terms"), instF))
        isEndo <- !(varNames %in% momNames)
        q <- length(momNames)
        tmpDat <- cbind(modelF, instF)
        add  <- survOptions$weights
        if (!is.null(vcovOptions$cluster))
            add <- cbind(as.matrix(vcovOptions$cluster), add)        
        if (!is.null(add))
            tmpDat <- cbind(tmpDat, add)
        na <- attr(get(na.action)(tmpDat), "na.action")[]
        if (!is.null(na))
        {
            modelF <- modelF[-na,,drop=FALSE]
            instF <- instF[-na,,drop=FALSE]
            if (!is.null(vcovOptions$cluster))
                {
                    if (is.null(dim(vcovOptions$cluster)))
                        vcovOptions$cluster <- vcovOptions$cluster[-na]
                    else
                        vcovOptions$cluster <- vcovOptions$cluster[-na,,drop=FALSE]
                }
            if (!is.null(survOptions$weights))
                survOptions$weights <- survOptions$weights[-na]
        }
        if (is.null(na))
            na <- integer()
        n <- nrow(modelF)
        list(modelF=modelF,  instF=instF, fRHS=fRHS, fLHS=fLHS, n=n, k=k, q=q,
             momNames=momNames, parNames=parNames, varNames=varNames, isEndo=isEndo,
             omit=na, vcovOptions=vcovOptions, survOptions=survOptions)
    }

.fModelData <- function(g, x, theta0, survOptions=list(), vcovOptions=list(),
                      na.action="na.omit", grad=NULL)
    {
        mom <- try(g(theta0, x))
        if (!is.null(grad))
        {
            dmom <- try(grad(theta0, x))
            if (inherits(dmom, "try-error"))
            {
                warning("grad could not be evaluated at the starting theta. Changed to NULL")
                grad <- NULL
            }
        }
        k <- length(theta0)        
        if (is.null(names(theta0)))
            {
                parNames <- paste("theta", 1:k, sep="")
                names(theta0) <- parNames
            } else {
                parNames <- names(theta0)
            }
        add  <- survOptions$weights
        if (!is.null(vcovOptions$cluster))
            add <- cbind(as.matrix(vcovOptions$cluster), add)        
        if (!is.null(add))
            {
                if (any(is.na(add)))
                    stop("weights or cluster contains missing values")
            }
        if (any(class(mom)=="try-error"))
            {
                msg <- paste("Cannot evaluate the moments at theta0\n",
                             attr(mom,"condition"))
                stop(msg)
            } else if (any(is.na(mom))) {
                stop("Some moments are NA's. Make sure you remove missing values from x")
            } else {
                q <-  ncol(mom)
                n <- nrow(mom)                
                if (!is.null(colnames(mom)))
                    momNames <- colnames(mom)
                else
                    momNames <- paste("h", 1:q, sep="")
            }
        list(q=q,n=n,k=k, momNames=momNames, parNames=parNames,
             varNames=character(), isEndo=logical(), omit=integer(),
             vcovOptions=vcovOptions, survOptions=survOptions, theta0=theta0,
             dfct=grad)
    }

.slModelData <- function(g,h,data, survOptions=list(), vcovOptions=list(),
                       na.action="na.omit")
    {
        res <- lapply(1:length(g), function(i) .lModelData(g[[i]], h[[i]], data,
                                                         list(), list(), "na.pass"))
        modelT <- lapply(res, function(x) attr(x$modelF, "terms"))
        instT <-  lapply(res, function(x) attr(x$instF, "terms"))
        allDat <-  do.call(cbind, lapply(res, function(x) cbind(x$modelF, x$instF)))
        allDat <- allDat[,!duplicated(colnames(allDat))]
        add  <- survOptions$weights
        if (!is.null(vcovOptions$cluster))
            add <- cbind(as.matrix(vcovOptions$cluster), add)        
        if (!is.null(add))
            tmpDat <- cbind(allDat, add)
        else
            tmpDat <- allDat
        na <- attr(get(na.action)(tmpDat), "na.action")[]
        if (!is.null(na))
        {
            allDat <- allDat[-na,,drop=FALSE]
            if (!is.null(vcovOptions$cluster))
                {
                    if (is.null(dim(vcovOptions$cluster)))
                        vcovOptions$cluster <- vcovOptions$cluster[-na]
                    else
                        vcovOptions$cluster <- vcovOptions$cluster[-na,,drop=FALSE]
                }
            if (!is.null(survOptions$weights))
                survOptions$weights <- survOptions$weights[-na]
        }
        if (is.null(na))
            na <- integer()
        parNames <- lapply(1:length(g), function(i) res[[i]]$parNames)
        momNames <- lapply(1:length(g), function(i) res[[i]]$momNames)
        isEndo <- lapply(1:length(g), function(i) res[[i]]$isEndo)
        varNames <- lapply(1:length(g), function(i) res[[i]]$varNames)
        k <- sapply(parNames, length)
        q <- sapply(momNames, length)
        n <- nrow(allDat)
        if (!is.null(names(g)))
            eqnNames=names(g)
        else
            eqnNames <- paste("Eqn", 1:length(g), sep="")
        list(data=allDat, modelT=modelT, instT=instT, parNames=parNames,
             momNames=momNames, k=k,q=q,n=n, eqnNames=eqnNames,
             varNames=varNames, isEndo=isEndo, omit=na,
             vcovOptions=vcovOptions, survOptions=survOptions)
    }

.snlModelData <- function(g,h,theta0, data, survOptions=list(), vcovOptions=list(),
                        na.action="na.omit")
    {
        res <- lapply(1:length(g), function(i) .nlModelData(g[[i]], h[[i]],
                                                          theta0[[i]], data, list(),
                                                          list(), "na.pass"))
        fRHS <- lapply(res, function(x) x$fRHS)
        fLHS <- lapply(res, function(x) x$fLHS)
        instT <-  lapply(res, function(x) attr(x$instF, "terms"))
        allDat <-  do.call(cbind, lapply(res, function(x) cbind(x$modelF, x$instF)))
        allDat <- allDat[,!duplicated(colnames(allDat))]
        add  <- survOptions$weights
        if (!is.null(vcovOptions$cluster))
            add <- cbind(as.matrix(vcovOptions$cluster), add)        
        if (!is.null(add))
            tmpDat <- cbind(allDat, add)
        else
            tmpDat <- allDat
        na <- attr(get(na.action)(tmpDat), "na.action")[]
        if (!is.null(na))
        {
            allDat <- allDat[-na,,drop=FALSE]
            if (!is.null(vcovOptions$cluster))
                {
                    if (is.null(dim(vcovOptions$cluster)))
                        vcovOptions$cluster <- vcovOptions$cluster[-na]
                    else
                        vcovOptions$cluster <- vcovOptions$cluster[-na,,drop=FALSE]
                }
            if (!is.null(survOptions$weights))
                survOptions$weights <- survOptions$weights[-na]
        }
        if (is.null(na))
            na <- integer()
        parNames <- lapply(1:length(g), function(i) res[[i]]$parNames)
        momNames <- lapply(1:length(g), function(i) res[[i]]$momNames)
        isEndo <- lapply(1:length(g), function(i) res[[i]]$isEndo)
        varNames <- lapply(1:length(g), function(i) res[[i]]$varNames)
        k <- sapply(parNames, length)
        q <- sapply(momNames, length)
        n <- nrow(allDat)
        if (!is.null(names(g)))
            eqnNames=names(g)
        else
            eqnNames <- paste("Eqn", 1:length(g), sep="")
        list(data=allDat, fRHS=fRHS, fLHS=fLHS, parNames=parNames,
             momNames=momNames, k=k,q=q,n=n, eqnNames=eqnNames, instT=instT,
             varNames=varNames, isEndo=isEndo, omit=na,
             vcovOptions=vcovOptions, survOptions=survOptions)
    }
