################### the main gel function ###################
########## The functions is to avoid having to build model objects

gel4 <- function (g, x=NULL, theta0=NULL,lambda0=NULL, getVcov=FALSE, 
                  gelType = c("EL","ET","EEL","HD", "REEL","ETEL","ETHD"),
                  vcov = c("MDS","iid","HAC"), grad=NULL,
                  vcovOptions=list(), centeredVcov = TRUE,
                  cstLHS=NULL, cstRHS=NULL, lamSlv=NULL,
                  rhoFct=NULL, initTheta=c("gmm", "theta0"),
                  data = parent.frame(),
                  coefSlv=c("optim","nlminb","constrOptim"),
                  smooth=FALSE,
                  lControl=list(), tControl=list())
{
    Call <- match.call()
    vcov <- match.arg(vcov)
    gelType <- match.arg(gelType)
    coefSlv <- match.arg(coefSlv)
    initTheta <- match.arg(initTheta)
    if (initTheta=="theta0" & is.null(theta0))
        stop("theta0 is required when initTheta='theta0'")
    model <- momentModel(g, x, theta0, grad, vcov, vcovOptions,
                         centeredVcov, data=data, smooth=smooth)
    if (initTheta == "theta0")
    {
        initTheta <- "modelTheta0"
        names(theta0) = model@parNames
    } else {
        theta0 <- NULL
    }
    if (!is.null(cstLHS))
    {
        model <- restModel(model, cstLHS, cstRHS)
        spec <- modelDims(model)
        if (!is.null(theta0))
            theta0 <- theta0[(names(theta0) %in% spec$parNames)]
    }
    fit <- gelFit(model=model, initTheta=initTheta, theta0=theta0,
                  lambda0=lambda0, vcov=getVcov, coefSlv=coefSlv,
                  lamSlv=lamSlv, tControl=tControl, lControl=lControl,
                  gelType=gelType, rhoFct=rhoFct)
    fit@call <- Call
    fit
}

