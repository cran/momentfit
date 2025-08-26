algoObj <- function(algo, start, fct, grad, solution, value, message, convergence)
{
    if (algo %in% c("optim", "nlminb", "constrOptim"))
    {
        obj <- switch(algo,
                      optim = new("minAlgoStd"),
                      nlminb = new("minAlgoStd",
                                   algo="nlminb", start="start", fct="objective",
                                   grad="gradient", value="objective"),
                      constrOptim = new("minAlgoStd", algo="constrOptim", start="theta",
                                        fct="f", grad="grad"))
        return(obj)
    } else if (algo == "nlm") {
        obj <- new("minAlgoNlm")
    } else {
        if (missing(start) | missing(fct))
            stop("You must provide the name of the arguments representing the function to minimize and the starting values")
        if (missing(value) | missing(solution) | missing(message) | missing(convergence))
            stop("You must provide the name of the output representing the solution, the function value, the convergence code and message.")
        obj <- new("minAlgoStd", algo=algo, start=start, fct=fct, grad=grad,
                   solution=solution,
                   value=value, message=message, convergence=convergence)
    }
    obj
}

setGeneric("minFit",
           def = function(object, start, fct, gr, ...) "Unknown algorithm")

setMethod("minFit", signature("minAlgoNlm"),
          function(object, start, fct, gr, ...)
          {
              solver <- object@algo
              arg <- list()
              if (missing(gr))
              {
                  f <- fct
              } else {
                  if (!is.function(gr))
                      stop("gr must be a function")
                  if (!isTRUE(all.equal(formals(fct), formals(gr))))
                      stop("Arguments in fct must be identical to arguments in gr")
                  f <- function()
                  {
                      arg <- as.list(match.call)[-1]
                      structure(do.call("fct", arg),
                                gradient=do.call("gr", arg))
                  }
                  formals(f) <- formals(fct)
              }
              arg[[object@fct]] <- fct
              arg[[object@start]] <- start
              arg <- c(arg, list(...))
              res <- do.call(solver, arg)
              ans <- list(solution = res[[object@solution]],
                          value = res[[object@value]])
              if (!is.na(object@convergence))
                  ans$convergence <- res[[object@convergence]]
              if (!is.na(object@message))
                  ans$message <- res[[object@message]]
              ans
          })



setMethod("minFit", signature("minAlgoStd"),
          function(object, start, fct, gr, ...)
          {
              solver <- object@algo
              arg <- list()
              arg[[object@fct]] <- fct
              if (!is.na(object@grad))
              {
                  if (!missing(gr))
                  {
                      if (!is.function(gr))
                          stop("gr must be a function")
                      if (!isTRUE(all.equal(formals(fct), formals(gr))))
                          stop("Arguments in fct must be identical to arguments in gr")
                      arg[[object@grad]] <- gr
                  }
              }
              arg[[object@start]] <- start
              arg <- c(arg, list(...))
              res <- do.call(solver, arg)
              ans <- list(solution = res[[object@solution]],
                          value = res[[object@value]])
              if (!is.na(object@convergence))
                  ans$convergence <- res[[object@convergence]]
              if (!is.na(object@message))
                  ans$message <- res[[object@message]]
              ans
          })



setMethod("print", "minAlgo",
          function(x, ...)
          {
              cat("Optimization algorithm\n")
              cat("**********************\n")
              cat("Name of the function: ", x@algo, "\n")
              invisible()
          })

setMethod("show", "minAlgo", function(object) print(object))
