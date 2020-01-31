######## Methods for summaryGmm object

## print
setMethod("print", "summaryGmm",
          function(x, digits=5, ...)
          {
              print(x@model)
              ntype <- matrix(c("Two-Step GMM", "Iterated GMM", "CUE",
                                "One-Step GMM with fixed weights","Two-Stage Least Squares",
                                "Evaluated at a fixed Theta; No estimation",
                                "One-Step Efficient M.D.E.",
                                "twostep","iter","cue","onestep","tsls", "eval","mde"),
                              ncol=2)              
              type <- ntype[match(x@type, ntype[, 2]), 1]
              spec <- modelDims(x@model)
              if (spec$q == spec$k) 
                  type <- "One-Step, Just-Identified"
              cat("\nEstimation: ", type, "\n")
              if (!is.null(x@convergence)) 
                  cat("Convergence Optim: ", x@convergence, "\n")
              if (!is.null(x@convIter)) 
                  cat("Convergence Iteration: ", x@convIter, "\n")
              if (x@type == "iter")
                  cat("Number of iterations: ", x@niter, "\n")
              if (length(x@wSpec) > 0)
                  {
                      if (is.numeric(x@model@vcovOptions$bw))
                          cat("Fixed Bandwidth: ", round(x@wSpec$bw, 3), "\n", sep="")
                      else
                          cat(x@model@vcovOptions$bw,
                              " Bandwidth: ", round(x@wSpec$bw, 3), "\n", sep="")
                  }
              if (x@breadOnly)
                  {
                      cat("vcov type: Bread \n")
                  } else {
                      cat("Sandwich vcov: ", x@sandwich, "\n", sep="")
                      if (x@sandwich && x@model@vcov == "MDS")
                          {
                              v <- ifelse(x@df.adj, "HC1", "HC0")
                              cat("Type of sandwich HCCM :", v, "\n", sep="")
                          } else if (x@sandwich && x@model@vcov == "HAC") {
                              cat("Type of sandwich HAC: as specified in the model definition\n")
                          }
                  }
              cat("coefficients:\n")
              printCoefmat(x@coef, digits=digits, ...)
              print(x@specTest)
              str <- x@strength
              if (!is.null(str$strength)) {
                  cat("\n", str$mess, "\n", sep = "")
                  str <- str$strength
                  for (i in 1:nrow(str)) cat(rownames(str)[i], ": F(", 
                                             str[i, 2], ", ", str[i, 3], ") = ", str[i, 1],
                                             " (P-Vavue = ", str[i, 4], ")\n")
              }})
## show
setMethod("show", "summaryGmm", function(object) print(object)) 

########## System of Equations

## print

setMethod("print", "summarySysGmm",
          function(x, digits = 5, ...) 
              {
                  print(x@model)
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
                                    "Seemingly Unrelated Regression", "twostep", "iter", 
                                    "cue", "onestep", "tsls", "eval", "EBEtwostep",
                                    "EBEiter", "EBEcue", "EBEonestep", "3SLS", "FIVE","SUR"),
                                  ncol = 2)
                  type <- ntype[match(x@type, ntype[, 2]), 1]
                  spec <- modelDims(x@model)
                  if (all(spec$q == spec$k)) 
                      type <- "One-Step, All equations are just-identified"
                  cat("\nEstimation: ", type, "\n")
                  EbE <- strtrim(x@type, 3) == "EBE"
                  if (!EbE)
                      {
                          if (!is.null(x@convIter)) 
                              cat("Convergence Iteration: ", x@convIter, "\n")
                          if (x@type == "iter") 
                              cat("Number of iterations: ", x@niter, "\n")
                          if (length(x@wSpec) > 0) {
                              if (is.numeric(x@model@vcovOptions$bw)) 
                                  cat("Fixed Bandwidth: ", round(x@wSpec$bw, 3), 
                                      "\n", sep = "")
                              else cat(x@model@vcovOptions$bw, " Bandwidth: ",
                                       round(x@wSpec$bw, 3), "\n", sep = "")
                              if (!is.null(x@convergence)) 
                                  cat("Convergence Optim: ", x@convergence, "\n")
                          }
                      }
                  if (x@breadOnly) {
                      cat("vcov type: Bread \n")
                  } else {
                      cat("Sandwich vcov: ", x@sandwich, "\n", sep = "")
                      if (x@sandwich && x@model@vcov == "MDS") {
                          v <- ifelse(x@df.adj, "HC1", "HC0")
                          cat("Type of sandwich HCCM :", v, "\n", sep = "")
                      }
                      else if (x@sandwich && x@model@vcov == "HAC") {
                          cat("Type of sandwich HAC: as specified in the model definition\n")
                      }
                  }
                  cat("coefficients:\n")
                  for (i in 1:length(x@coef))
                      {
                          sleg <- i==length(x@coef)
                          str <- x@strength[[i]]
                          cat("\n",names(x@coef)[i], ": \n", sep="")
                          if (EbE)
                              {
                                  if (!is.null(x@convIter)) 
                                      cat("Convergence Iteration: ", x@convIter[i], "\n")
                                  if (x@type == "EBEiter") 
                                      cat("Number of iterations: ", x@niter[i], "\n")
                                  if (length(x@wSpec) > 0) {
                                      if (is.numeric(x@model@bw)) 
                                          cat("Fixed Bandwidth: ", round(x@wSpec$bw, 3), 
                                              "\n", sep = "")
                                      else cat(x@model@bw, " Bandwidth: ",
                                               round(x@wSpec$bw, 3), "\n", sep = "")
                                  }
                                  if (!is.null(x@convergence)) 
                                      cat("Convergence Optim: ", x@convergence[i], "\n")
                              }
                          printCoefmat(x@coef[[i]], digits = digits, signif.legend=sleg, ...)
                          if (!is.null(str$strength)) {
                              cat("\n", str$mess, "\n", sep = "")
                              str <- str$strength
                              for (i in 1:nrow(str))
                                  cat(rownames(str)[i], ": F(", str[i, 2], ", ",
                                      str[i, 3], ") = ", str[i, 1],  " (P-Vavue = ",
                                      str[i, 4], ")\n")
                          }
                      }                          
                  print(x@specTest)
              })

## show

setMethod("show", "summarySysGmm", function(object) print(object))


