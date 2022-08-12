################### the main gmm functions ###################
########## These functions ar to avoid having to builf model objects

gmm4 <- function (g, x, theta0 = NULL, grad = NULL, 
                  type = c("twostep", "iter", "cue", "onestep"),
                  vcov = c("iid", "HAC", "MDS", "TrueFixed", "CL"),
                  initW = c("ident", "tsls", "EbyE"), weights = "optimal", 
                  itermaxit = 50, cstLHS=NULL, cstRHS=NULL,
                  vcovOptions=list(),survOptions=list(),
                  itertol = 1e-07, centeredVcov = TRUE,
                  data = parent.frame(), ...) 
{
    Call <- match.call()
    vcov <- match.arg(vcov)
    type <- match.arg(type)
    initW <- match.arg(initW)
    if (vcov == "TrueFixed")
    {
        if (!is.matrix(weights) ||
            !inherits(weights,c("gmmWeights", "sysGmmWeigths")))
            stop("With TrueFixed vcov the weights must be provided")
        efficientWeights <- TRUE
        vcov2 <- "iid"
    } else {
        efficientWeights <- FALSE
        vcov2 <- vcov
    }
    if (is.list(g))
    {
        ## Formula or sysGMM? Need to find a better way.
        model <- NULL
        if (is.null(x) & !is.null(theta0))
            model <- try(momentModel(g=g, x=x, theta0=theta0, grad=grad, vcov=vcov2,
                                  vcovOptions=vcovOptions,survOptions=survOptions,
                                  centeredVcov=centeredVcov, data=data), silent=TRUE)
        if (is.null(model) || inherits(model,"try-error"))
            model <- sysMomentModel(g=g, h=x, theta0=theta0, vcov=vcov2,
                                 vcovOptions=vcovOptions,survOptions=survOptions,
                                 centeredVcov=centeredVcov, data=data)
    } else {
        model <- momentModel(g=g, x=x, theta0=theta0, grad=grad, vcov=vcov2,
                          vcovOptions=vcovOptions,survOptions=survOptions,
                          centeredVcov=centeredVcov, data=data)
        if (initW == "EbyE")
        {
            warning("initW cannot be EbyE for single equations, initW set to ident")
            initW="ident"
        }
    }
    if (!is.null(cstLHS))
        model <- restModel(model, cstLHS, cstRHS)
    fit <- gmmFit(model=model, type=type, itertol=itertol, initW=initW,
                  weights=weights, itermaxit=itermaxit,
                  efficientWeights=efficientWeights, ...)
    fit@call <- Call
    fit
}


setMethod("tsls", "formula",
          function(model, x, vcov = c("iid", "HAC", "MDS", "CL"),                   
                   vcovOptions=list(), survOptions=list(), centeredVcov = TRUE,
                   data = parent.frame())
          {
              Call <- match.call(call=sys.call(sys.parent()-1L))
              vcov <- match.arg(vcov)
              model <- momentModel(g = model, x = x, vcov = vcov,
                                vcovOptions=vcovOptions,survOptions=survOptions,
                                centeredVcov = centeredVcov, data = data)
              obj <- tsls(model)
              obj@call <- Call
              obj
              })

setMethod("tsls", "list",
          function(model, x=NULL, vcov = c("iid", "HAC", "MDS", "CL"),
                   vcovOptions=list(), survOptions=list(), centeredVcov = TRUE,
                   data = parent.frame())
          {
              Call <- match.call(call=sys.call(sys.parent()-1L))              
              vcov <- match.arg(vcov)
              model <- sysMomentModel(g = model, h = x, vcov = vcov,
                                   vcovOptions=vcovOptions,survOptions=survOptions,
                                   centeredVcov = centeredVcov, data = data)
              obj <- tsls(model)
              obj@call <- Call
              obj
              })


setMethod("ThreeSLS", "list",
          function(model, x=NULL, vcov = c("iid", "HAC", "MDS", "CL"),
                   vcovOptions=list(), survOptions=list(), centeredVcov = TRUE,
                   data = parent.frame())
          {
              Call <- match.call(call=sys.call(sys.parent()-1L))              
              vcov <- match.arg(vcov)
              model <- sysMomentModel(g = model, h = x, vcov = vcov,
                                   vcovOptions=vcovOptions,survOptions=survOptions,
                                   centeredVcov = centeredVcov, data = data)
              obj <- ThreeSLS(model)
              obj@call <- Call
              obj
              })
