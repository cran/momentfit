
##################  Constructor for the sysGmmModels classes   #####################

.sysFct <- function(g, h=NULL, theta0=NULL,grad=NULL,
                    vcov = c("iid", "HAC", "MDS", "CL"),
                    vcovOptions=list(), centeredVcov = TRUE, data=parent.frame(),
                    na.action="na.omit", survOptions=list())
{
    vcov <- match.arg(vcov)
    if (vcov == "iid")
        vcov <- "MDS"
    if (!is.list(g))
        stop("g must be a list")
    if (!all(sapply(g, function(gi)
        inherits(gi, "function"))))
        stop("g must be a list of functions")
    if (!is.list(theta0))
        stop("theta0 must be a list")
    if (length(g) != length(theta0))
        stop("The length of g and theta0 must be the same")
    neqn <- length(g)
    if (is.null(grad))
    {
        grad <- lapply(1:neqn, function(i) NULL)
    } else {
        if (!is.list(grad))
            stop("grad must be a list")
        if (length(grad) != neqn)
            stop("The length is grad must be equal to the length of g")
    }
    if (is.null(names(g)))
    {
        eqnNames <- paste("Eqn", 1:neqn, sep="")
    } else {
        eqnNames <- names(g)
    }
    mod <- lapply(1:neqn, function(i)
        momentModel(g=g[[i]], x=h, theta0=theta0[[i]],grad=grad[[i]],
                    vcov=vcov, vcovOptions=vcovOptions,
                    centeredVcov=centeredVcov, data=data,
                    na.action=na.action, survOptions=survOptions,
                    smooth=FALSE))
    new("sfunctionModel",
        X=mod[[1]]@X,
        fct=lapply(mod, function(mi) mi@fct),
        dfct=lapply(mod, function(mi) mi@dfct),        
        vcov=mod[[1]]@vcov,
        theta0=lapply(mod, function(mi) mi@theta0),
        n=mod[[1]]@n,
        q=sapply(mod, function(mi) mi@q),
        k=sapply(mod, function(mi) mi@k),
        parNames=lapply(mod, function(mi) mi@parNames),
        momNames=lapply(mod, function(mi) mi@momNames),
        eqnNames=eqnNames, vcovOptions=mod[[1]]@vcovOptions,
        centeredVcov=mod[[1]]@centeredVcov,
        sameMom=FALSE, SUR=FALSE,
        varNames=lapply(mod, function(mi) mi@varNames), 
        omit=mod[[1]]@omit, survOptions=mod[[1]]@survOptions,
        sSpec=mod[[1]]@sSpec, smooth=mod[[1]]@smooth)
}


sysMomentModel <- function(g, h=NULL, theta0=NULL,grad=NULL,
                           vcov = c("iid", "HAC", "MDS", "CL"),
                           vcovOptions=list(), centeredVcov = TRUE,
                           data=parent.frame(),
                           na.action="na.omit", survOptions=list())
{
    vcov <- match.arg(vcov)
    if (!is.list(vcovOptions) | !is.list(survOptions))
        stop("vcovOptions and survOptions must be a list")
    vcovOptions <- .getVcovOptions(vcov, data, vcovOptions, FALSE)
    survOptions <- .getSurvOptions(data, survOptions)
    if (!is.list(data) && !is.environment(data)) 
        stop("'data' must be a list or an environment")
    if (!is.list(g))
        stop("For system of equations, g must be lists of formulas")
    if (all(sapply(g, function(gi) inherits(gi, "function"))))
    {
        mod <- .sysFct(g, h, theta0,grad, vcov,
                    vcovOptions, centeredVcov, data,
                    na.action="na.omit", survOptions=list())
        return(mod)
    }
    clg <- sapply(g, class)
    if (!all(clg=="formula"))
        stop("For system of equations, g must be lists of formulas")
    if (length(g) == 1)
        stop("There is only one equation. Use momentModel() instead")
    ## Try to see if parameters are in formulas
    if (!is.null(theta0))
    {
        if (!is.list(theta0))
            stop("For system of equations, theta0 must be a list of named vectors")
        if (length(theta0) != length(g))
            stop("The length of theta0 must match the length of g and h")
        chk <- sapply(1:length(theta0),
                      function(i) isTRUE(all(names(theta0[[i]]) %in% all.vars(g[[i]]))))
        if (all(chk))
        {
            nonlin <- TRUE
            chk2 <- sapply(theta0[-1],
                           function(l) any(names(l) %in% names(theta0[[1]])))
            if (any(chk2))
                stop("Coefficient names across equations must be different")
        } else if (all(!chk)){
            nonlin <- FALSE
        } else {
            stop("g must all be linear or all nonlinear")
        }
    } else {
        nonlin <- FALSE
    }
    varg <- lapply(1:length(g), function(i) .formSpec(g[[i]], theta0[[i]]))
    if (is.null(h))
    {
        SUR <- TRUE
        sameMom <- TRUE
        anyEndo <- rep(FALSE, length(g))
        chk <- sapply(varg, function(v) length(v$lhs)==0)
        if (any(chk))
            stop("For SUR, there must be a LHS")
        reg <- lapply(1:length(g), function(i) varg[[i]]$rhs)
        reg <- unique(do.call("c", reg)) 
        if (!nonlin)
            intercept <- any(sapply(varg, function(v) v$intercept == 1))
        else
            intercept <- TRUE
        reg <- paste(reg, collapse="+", sep="")
        reg <- paste("~", reg, ifelse(intercept, "", "-1"), sep="")
        reg <- as.formula(reg, .GlobalEnv)
        h <- rep(list(reg), length(g))
    } else {
        SUR <- FALSE
        if (inherits(h,"formula"))
        {
            h <- rep(list(h), length(g))
            sameMom <- TRUE
        } else if (!is.list(h)) {
            stop("h must be a list or a formula")
        } else if (length(h) == 1) {
            h <- rep(h, length(g))
            sameMom <- TRUE
        } else if (length(h) != length(g)) {
            stop("With different instruments, we must have length(h)= length(g)")
        } else {
            sameMom <- FALSE
        }
        varh <- lapply(h, .formSpec)
        anyEndo <- !sapply(1:length(g),
                           function(i) all(varg[[i]]$rhs %in% varh[[i]]$rhs) &
                                       (varg[[i]]$intercept==varh[[i]]$intercept))
    }
    if (!nonlin)
    {
        model <- .slModelData(g,h,data, survOptions, vcovOptions, na.action)
        isEndo <- lapply(1:length(g),
                         function(i) model$parNames[[i]] %in% model$momNames[[i]])
        gmodel <- new("slinearModel", data=model$data, 
                      instT=model$instT, modelT=model$modelT,
                      vcov=vcov, vcovOptions=model$vcovOptions,
                      centeredVcov = centeredVcov, k=model$k,
                      q=model$q, n=model$n, parNames=model$parNames,
                      momNames=model$momNames, eqnNames=model$eqnNames,
                      sameMom=sameMom, SUR=SUR, varNames=model$varNames,
                      isEndo=model$isEndo, omit=model$omit,
                      survOptions=model$survOptions, smooth=FALSE)
    } else {
        model <- .snlModelData(g, h, theta0, data, survOptions, vcovOptions, na.action)
        gmodel <- new("snonlinearModel", data=model$data, instT=model$instT,
                      theta0=theta0,fRHS=model$fRHS,eqnNames=model$eqnNames,
                      fLHS=model$fLHS, vcov=vcov, vcovOptions=model$vcovOptions,
                      centeredVcov = centeredVcov, k=model$k, q=model$q,
                      n=model$n, parNames=model$parNames,
                      momNames=model$momNames, sameMom=sameMom, SUR=SUR,
                      varNames=model$varNames, isEndo=model$isEndo, omit=model$omit,
                      survOptions=model$survOptions, smooth=FALSE)
    }
    gmodel
}


.formSpec <- function(f, theta=NULL)
{
    if (length(f) == 3)
    {  
        rhs <- all.vars(f[[3]])
        if (!is.null(theta))
            rhs <- rhs[!(rhs%in%names(theta))]
        lhs <- all.vars(f[[2]])
        if (!is.null(theta))
            lhs <- lhs[!(lhs%in%names(theta))]
    } else {
        lhs <- character()
        rhs <- all.vars(f[[2]])
        if (!is.null(theta))
            rhs <- rhs[!(rhs%in%names(theta))]
    }
    list(rhs=rhs, lhs=lhs, intercept=attr(terms(f), "intercept"))
}
