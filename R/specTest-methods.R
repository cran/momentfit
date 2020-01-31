### print

setMethod("print", "specTest", 
          function (x, digits = 5, ...) 
          {
              cat("\n", x@testname, "\n")
              print.default(x@test, digits = digits, print.gap = 2, 
                            quote = FALSE)
              cat("\n")
          })

### show

setMethod("show", "specTest", function(object) print(object)) 
