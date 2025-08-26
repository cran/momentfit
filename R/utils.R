MPinv <- function(x)
{
    d <- dim(x)
    tr <- FALSE
    if (d[1]<d[2])
    {
        x <- t(x)
        d <- d[c(2,1)]
        tr <- TRUE
    }
    res <- .Fortran(F_pinv, as.double(x), as.integer(d[1]), as.integer(d[2]),
                    r=integer(1), L=double(d[2]^2))
    L <- matrix(res$L, d[2], d[2])[,1:res$r]
    if (dim(L)[1] != dim(L)[2])
    {
        M <- solve(crossprod(L))
        ix <- crossprod(M%*%t(L))%*%t(x)
    } else {
        ix <- crossprod(forwardsolve(L, diag(dim(L)[1])))%*%t(x)
    }
    if (tr)
        ix = t(ix)
    ix
}
