## print and show

setMethod("print", "hypothesisTest",
          function(x, ...) {
              cat(x@type, "\n")
              cat("***********\n")
              cat("The Null Hypothesis:\n")
              for (i in 1:length(x@hypothesis)) cat("\t", x@hypothesis[i], "\n")
              cat("Distribution: ", x@dist, " with ", x@df, " degrees of freedom\n", sep="")
              print(data.frame(Statistics=x@test, Pvalue=x@pvalue), ...)
          })

setMethod("show", "hypothesisTest", function(object) print(object))

#########
