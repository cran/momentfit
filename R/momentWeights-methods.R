####### Methods for gmmWeights objects
#####################################


################ method to compute X'WX where W is the object weights #############

setGeneric("quadra", function(w, x, y, ...) standardGeneric("quadra"))

setMethod("quadra", c("momentWeights", "matrixORnumeric", "missing"),
          function(w, x, y) {
              x <- as.matrix(x)
              if (is.character(w@w))
                  {
                      obj <- crossprod(x)
                  } else {
                      if (w@type == "weights")
                          {
                              obj <- crossprod(x,w@w)%*%x
                          } else if (w@type == "vcov") {
                              obj <- crossprod(x, solve(w@w,x))
                          } else if (w@type == "chol") {
                              obj <- crossprod(forwardsolve(t(w@w),x))
                          } else {
                              R <- qr.R(w@w)
                              if (!all(w@w$pivot==(1:ncol(R))))
                                  {
                                      iw <- matrix(NA, ncol(R), ncol(R))
                                      iw[w@w$pivot, w@w$pivot] <- chol2inv(R)
                                      obj <- crossprod(x,iw)%*%x
                                  } else {
                                      obj <- crossprod(forwardsolve(t(R),x))
                                  }
                          }
                  }
              if (all(dim(obj)==1))
                  c(obj)
              else
                  obj                  
          })

setMethod("quadra", c("momentWeights", "matrixORnumeric", "matrixORnumeric"),
          function(w, x, y) {
              x <- as.matrix(x)
              y <- as.matrix(y)
              if (is.character(w@w))
              {                  
                  obj <- crossprod(x,y)
              } else {
                  if (w@type == "weights")
                      {
                          obj <- crossprod(x,w@w)%*%y
                      } else if (w@type == "vcov") {
                          obj <- crossprod(x, solve(w@w,y))
                      } else if (w@type == "chol") {
                          T1 <- forwardsolve(t(w@w), x)
                          T2 <- forwardsolve(t(w@w), y)
                          obj <- crossprod(T1,T2)
                      } else {
                          R <- qr.R(w@w)
                          if (!all(w@w$pivot==(1:ncol(R))))
                              {
                                  iw <- matrix(NA, ncol(R), ncol(R))
                                  iw[w@w$pivot, w@w$pivot] <- chol2inv(R)
                                  obj <- crossprod(x,iw)%*%y
                              } else {
                                  T1 <- forwardsolve(t(R), x)
                                  T2 <- forwardsolve(t(R), y)
                                  obj <- crossprod(T1,T2)
                              }
                      }
              }
              if (all(dim(obj)==1))
                  c(obj)
              else
                  obj                  
          })

setMethod("quadra", c("momentWeights", "missing", "missing"),
          function(w, x, y) {
              if (is.character(w@w))
              {                  
                  obj <- "Identity"
              } else {
                  if (w@type == "weights")
                      {
                          obj <- w@w
                      } else if (w@type == "vcov") {
                          obj <- solve(w@w)
                      } else if (w@type == "chol") {
                          obj <- chol2inv(w@w)
                      } else {
                          R <- qr.R(w@w)
                          if (!all(w@w$pivot==(1:ncol(R))))
                              {
                                  iw <- matrix(NA, ncol(R), ncol(R))
                                  iw[w@w$pivot, w@w$pivot] <- chol2inv(R)
                                  obj <- iw
                              } else {
                                  obj <- chol2inv(R)
                              }
                      }
              }
              if (all(dim(obj)==1))
                  c(obj)
              else
                  obj                  
          })


setMethod("print", "momentWeights",
          function(x, ...){
              cat("Moment weights matrix object\n")
              print(quadra(x), ...)})

setMethod("show", "momentWeights", function(object) print(object))

### Subsetting method '['


setMethod("[", signature("momentWeights","numeric","missing"),
          function(x, i, j)
              {
                  i <- as.integer(i)
                  if (is.character(x@w))
                      return(x)
                  w <- quadra(x)[i,i]                 
                  x@w <- w
                  x@type <- "weights"
                  x
              })

setMethod("[", signature("momentWeights","missing","missing"),
          function(x, i,j)
              {
                  if (is.character(x@w))
                      return(x)
                  w <- quadra(x)                 
                  x@w <- w
                  x@type <- "weights"
                  x
              })


############## The following is for system weights
###################################################

### quadra

setMethod("quadra", c("sysMomentWeights", "matrixORnumeric", "missing"),
          function(w, x, y) {
              x <- as.matrix(x)
              if (is.character(w@w))
                  {
                      obj <- crossprod(x)
                  } else {
                      q <- sapply(w@momNames, length)
                      if (w@type == "weights")
                          {
                              obj <- crossprod(x,w@w)%*%x
                          } else if (w@type == "vcov") {
                              obj <- crossprod(x, solve(w@w,x))
                          } else if (w@type == "iid") {
                              if (w@sameMom)
                              {
                                  Rw <- qr.R(w@w)                                 
                                  if (!all(w@w$pivot==(1:ncol(Rw))))
                                  {
                                      iw <- matrix(NA, ncol(Rw), ncol(Rw))
                                      iw[w@w$pivot, w@w$pivot] <- chol2inv(Rw)
                                  } else {
                                      iw <- chol2inv(Rw)
                                  }
                                  isig <- chol2inv(w@Sigma)
                                  v <- kronecker(isig, iw)
                                  obj <- crossprod(x,v)%*%x
                              } else {
                                  q <- sapply(w@momNames, length)
                                  v <- .SigmaZZ(w@w, w@Sigma, q)
                                  obj <- crossprod(x, solve(v,x))
                              }
                          } else {
                              R <- qr.R(w@w)
                              if (!all(w@w$pivot==(1:ncol(R))))
                                  {
                                      iw <- matrix(NA, ncol(R), ncol(R))
                                      iw[w@w$pivot, w@w$pivot] <- chol2inv(R)
                                      obj <- crossprod(x,iw)%*%x
                                  } else {
                                      obj <- crossprod(forwardsolve(t(R),x))
                                  }
                          }
                  }
              if (all(dim(obj)==1))
                  c(obj)
              else
                  obj                  
          })

setMethod("quadra", c("sysMomentWeights", "matrixORnumeric", "matrixORnumeric"),
          function(w, x, y) {
              x <- as.matrix(x)
              y <- as.matrix(y)
              if (is.character(w@w))
                  {
                      obj <- crossprod(x,y)
                  } else {
                      q <- sapply(w@momNames, length)
                      if (w@type == "weights")
                          {
                              obj <- crossprod(x,w@w)%*%y
                          } else if (w@type == "vcov") {
                              obj <- crossprod(x, solve(w@w,y))
                          } else if (w@type == "iid") {
                              if (w@sameMom)
                              {
                                  Rw <- qr.R(w@w)                                 
                                  if (!all(w@w$pivot==(1:ncol(Rw))))
                                  {
                                      iw <- matrix(NA, ncol(Rw), ncol(Rw))
                                      iw[w@w$pivot, w@w$pivot] <- chol2inv(Rw)
                                  } else {
                                      iw <- chol2inv(Rw)
                                  }
                                  isig <- chol2inv(w@Sigma)
                                  v <- kronecker(isig, iw)
                                  obj <- crossprod(x,v)%*%y
                              } else {
                                  q <- sapply(w@momNames, length)
                                  v <- .SigmaZZ(w@w, w@Sigma, q)
                                  obj <- crossprod(x, solve(v,y))
                              }
                          } else {
                              R <- qr.R(w@w)
                              if (!all(w@w$pivot==(1:ncol(R))))
                                  {
                                      iw <- matrix(NA, ncol(R), ncol(R))
                                      iw[w@w$pivot, w@w$pivot] <- chol2inv(R)
                                      obj <- crossprod(x,iw)%*%y
                                  } else {
                                      T1 <- forwardsolve(t(R), x)
                                      T2 <- forwardsolve(t(R), y)
                                      obj <- crossprod(T1,T2)
                                  }
                          }
                  }
              if (all(dim(obj)==1))
                  c(obj)
              else
                  obj                  
          })

setMethod("quadra", c("sysMomentWeights", "missing", "missing"),
          function(w, x, y) {
              if (is.character(w@w))
              {                  
                  obj <- "Identity"
              } else {
                  if (w@type == "weights")
                      {
                          obj <- w@w
                      } else if (w@type == "vcov") {
                          obj <- solve(w@w)
                      } else if (w@type == "iid") {
                          if (w@sameMom)
                          {
                              Rw <- qr.R(w@w)                                 
                              if (!all(w@w$pivot==(1:ncol(Rw))))
                              {
                                  iw <- matrix(NA, ncol(Rw), ncol(Rw))
                                  iw[w@w$pivot, w@w$pivot] <- chol2inv(Rw)
                              } else {
                                  iw <- chol2inv(Rw)
                              }
                              isig <- chol2inv(w@Sigma)
                              obj <- kronecker(isig, iw)
                          } else {
                                  q <- sapply(w@momNames, length)
                                  v <- .SigmaZZ(w@w, w@Sigma, q)
                                  obj <- solve(v)
                          }
                      } else {
                          R <- qr.R(w@w)
                          if (!all(w@w$pivot==(1:ncol(R))))
                              {
                                  iw <- matrix(NA, ncol(R), ncol(R))
                                  iw[w@w$pivot, w@w$pivot] <- chol2inv(R)
                                  obj <- iw
                              } else {
                                  obj <- chol2inv(R)
                              }
                      }
              }
              if (all(dim(obj)==1))
                  c(obj)
              else
                  obj                  
          })


setMethod("print", "sysMomentWeights",
          function(x, ...){
              cat("Moment weights matrix object\n")
              w <- quadra(x)
              if (is.matrix(w))
              {
                  q <- sapply(x@momNames, length)
                  wn <- paste(rep(x@eqnNames, q), ".", do.call("c", x@momNames), sep="")
                  dimnames(w) <-list(wn,wn)
              }
              print(w, ...)})

setMethod("show", "sysMomentWeights", function(object) print(object))

### Subsetting method '['

setMethod("[", signature("sysMomentWeights","numeric", "missing"),
          function(x, i, j)
              {
                  i <- as.integer(i)
                  if (is.character(x@w))
                      return(x)
                  q <- sapply(x@momNames, length)
                  eq <- rep(1:length(x@eqnNames), q)
                  weq <- eq %in% i
                  w <- quadra(x)[weq,weq]
                  eqnNames <- x@eqnNames[(1:length(x@eqnNames) %in% i)]
                  x@type <- "weights"
                  x@w <- w
                  x@momNames <- x@momNames[i]
                  x@Sigma <- NULL
                  x@eqnNames <- eqnNames
                  x
              })

setMethod("[", signature("sysMomentWeights","missing", "list"),
          function(x, i, j)
          {
              if (is.character(x@w))
                  return(x)
              q <- sapply(x@momNames, length)
              if (length(j) != length(q))
                  stop("j must be a list with a length equals to the number of equations")
              sel <- vector()
              for (l in 1:length(j))
              {
                      if (length(j[[l]]) > 0)
                      {                         
                              q2 <- q[l]
                              if (!all(abs(j[[l]]) %in% (1:q2))) 
                                  stop("SubMoment must be between 1 and q")
                              x@momNames[[l]] <- x@momNames[[l]][j[[l]]]
                              sel <- c(sel, 1:q2 %in% j[[l]])
                      } else {
                          sel <- c(sel, rep(TRUE, q2))
                      }
              }
              w <- quadra(x)[sel,sel]
              x@type <- "weights"
              x@w <- w
              x@Sigma <- NULL
              x
          })

setMethod("[", signature("sysMomentWeights", "numeric", "list"),
          function(x, i, j) x[i][,j])
