#' Within transformation with unbalanced data
#'
#' Either \code{i}, \code{j}, \code{t} and \code{value} parameters or \code{m} data frame with all these four vectors are required.
#' @param i country 1 (vector)
#' @param j country 2 (vector)
#' @param t year (vector)
#' @param value vector to be transformed
#' @param m a \code{matrix} holding 4 columns for \code{i}, \code{j}, \code{t} and \code{value}
#' @return vector
#' @references László Balázsi, László Mátyás and Tom Wansbeek (2014): The Estimation of Multi-dimensional Fixed Effects Panel Data Models. Formula (22) and (23).
#' @export
#' @importFrom Matrix t rowSums colSums
#' @aliases WithinTransformation1c WithinTransformation2c WithinTransformation3c WithinTransformation4c WithinTransformation5c WithinTransformation6c
WithinTransformation1c <- function(i = m[, 1], j = m[, 2], t = m[, 3], value = m[, 4:ncol(m)], m)  { # formula (22)

    M <- cbind(i, j, t, value)

    ## TODO: check if i and j are factors (in other functions as well)
    n  <- length(levels(i))
    nt <- length(unique(t))

    v <- sparseMatrix(1:n^2, 1:n^2, x = 1)
    a <- apply(expand.grid(levels(i), levels(i)), 1, paste, collapse = '-')
    Vs <- lapply(setNames(unique(t), unique(t)), function(x) {
        w <- which(t == x)
        r <- which(a %in% paste(i[w], j[w], sep = '-'))
        structure(v[r, ], .Names = x)
    })

    D1 <- do.call(rBind, Vs)

    l <- sapply(Vs, function(x) x@Dim[1])
    names(l) <- 1:length(l)
    js <- as.numeric(rep(names(l), times = l))

    D2 <- sparseMatrix(1:length(js), js, x = 1)

    dN2 <- Matrix::t(D1) %*% D1 ## Matrix::tcrossprod

    ## inverse of dN2 with reduced number of cells
    ix <- union(which(Matrix::rowSums(dN2) > 0), which(Matrix::colSums(dN2) > 0))
    dN2i <- dN2
    dN2i[ix, ix] <- solve(dN2[ix, ix])

    Aa  <- Matrix::t(D2) %*% D1
    Da <- D2 - (D1 %*% dN2i %*% Matrix::t(Aa))


    dTx <- Matrix::t(D2) %*% D2
    It  <- sparseMatrix(1:dim(D1)[1], 1:dim(D1)[1], x = 1)
    Qa <- Matrix::t(D2) %*% Da

    ## Moore–Penrose pseudoinverse of Qa
    if (require('rfunctions')) {
        Qai <- geninv(as.matrix(Qa))
    } else {
        Qai <- ginv(as.matrix(Qa))
    }
    Qai <- Matrix(Qai)

    ## this should work on smaller matrices (formula 22)
    ## Pa <- (It - (D1 %*% dN2i %*% t(D1))) - (Da %*% Qai %*% t(Da))

    fia <- Qai %*% Matrix::t(Da) %*% value

    t <- as.numeric(factor(t))
    i <- as.numeric(i)
    j <- as.numeric(j)

    ## res <- simplify2array(mclapply(1:length(value), function(index) {
    ##     value[index] - mean(value[which(i == i[index] & j == j[index])]) - fia[as.numeric(t[index])]
    ## }))
    res <- sapply(1:length(value), function(index) {
        value[index] - mean(value[which(i == i[index] & j == j[index])]) - fia[as.numeric(t[index])]
    })

    res

}
#' @export
WithinTransformation2c <- WithinTransformation2
#' @export
WithinTransformation3c <- WithinTransformation3
#' @export
WithinTransformation4c <- WithinTransformation4
#' @export
WithinTransformation5c <- WithinTransformation5
#' @export
WithinTransformation6c <- function(i = m[, 1], j = m[, 2], t = m[, 3], value = m[, 4:ncol(m)], m)  { # formula (23)

    M <- cbind(i, j, t, value)

    n  <- length(levels(i))
    nt <- length(unique(t))

    v <- sparseMatrix(1:n^2, 1:n^2, x = 1)
    a <- apply(expand.grid(levels(i), levels(i)), 1, paste, collapse = '-')
    Vs <- lapply(setNames(unique(t), unique(t)), function(x) {
        w <- which(t == x)
        r <- which(a %in% paste(i[w], j[w], sep = '-'))
        structure(v[r, ], .Names = x)
    })
    D1 <- do.call(rBind, Vs)

    In <- sparseMatrix(1:n, 1:n, x = 1)
    ln <- rep(1, n)
    Us <- lapply(unique(t), function(x) {
        w <- which(t == x)
        r <- which(a %in% paste(i[w], j[w], sep = '-'))
        kronecker(In, ln)[r, ]
    })
    Ws <- lapply(unique(t), function(x) {
        w <- which(t == x)
        r <- which(a %in% paste(i[w], j[w], sep = '-'))
        kronecker(ln, In)[r, ]
    })

    l  <- sapply(Us, function(x) x@Dim)
    l  <- cbind(c(1,1), l)
    D3 <- D2 <- Matrix(0, nrow = sum(l[1, ])-1, ncol = sum(l[2, ])-1)
    for (index in 1:nt) {
        rstart <- cumsum(l[1, ])[index]
        cstart <- cumsum(l[2, ])[index]
        D2[rstart:(rstart-1+l[1, ][index+1]), cstart:(cstart-1+l[2, ][index+1])] <- Us[[index]]
        D3[rstart:(rstart-1+l[1, ][index+1]), cstart:(cstart-1+l[2, ][index+1])] <- Ws[[index]]
    }

    ## helper matrix for inverse
    D1tD1 <- Matrix::t(D1) %*% D1
    ix <- union(which(Matrix::rowSums(D1tD1) > 0), which(Matrix::colSums(D1tD1) > 0))
    D1tD1i <- D1tD1
    D1tD1i[ix, ix] <- solve(D1tD1[ix, ix])

    It <- sparseMatrix(1:dim(D1)[1], 1:dim(D1)[1], x = 1)
    B <- It - (D1 %*% D1tD1i %*% Matrix::t(D1))

    ## another helper matrix
    BD2ti <- BD2t <- Matrix::t(B %*% D2) %*% (B %*% D2)
    ix <- union(which(Matrix::rowSums(BD2t) > 0), which(Matrix::colSums(BD2t) > 0))
    if (require('rfunctions')) {
        BD2ti[ix, ix] <- geninv(as.matrix(BD2ti[ix, ix]))
    } else {
        BD2ti[ix, ix] <- ginv(as.matrix(BD2ti[ix, ix]))
    }

    C <- B - (B %*% D2) %*% BD2ti %*% Matrix::t(B %*% D2)

    Db <- (It - D1 %*% D1tD1i %*% Matrix::t(D1)) %*% D2
    Qb <- Matrix::t(D2) %*% Db

    ## Moore–Penrose pseudoinverse of Qb
    ix <- union(which(Matrix::rowSums(Qb) > 0), which(Matrix::colSums(Qb) > 0))
    Qbi <- Qb
    if (require('rfunctions')) {
        Qbi[ix, ix] <- geninv(as.matrix(Qb[ix, ix]))
    } else {
        Qbi[ix, ix] <- ginv(as.matrix(Qb[ix, ix]))
    }

    fib <- Qbi %*% Matrix::t(Db) %*% value

    Qi <- Q <- Matrix::t(C %*% D3) %*% (C %*% D3)

    ## Qb and the Moore–Penrose pseudoinverse of Qb
    ix <- union(which(Matrix::rowSums(Q) > 0), which(Matrix::colSums(Q) > 0))
    if (require('rfunctions')) {
        Qi[ix, ix] <- geninv(as.matrix(Q[ix, ix]))
    } else {
        Qi[ix, ix] <- ginv(as.matrix(Q[ix, ix]))
    }

    ## omega
    o <- Qi %*% Matrix::t(C %*% D3) %*% value

    ## xi
    k <- Qbi %*% Matrix::t(Db) %*% D3 %*% o

    Ab <- Matrix::t(D2) %*% D1
    A  <- Matrix::t(D3) %*% D1

    aij <- apply(expand.grid(seq_len(n), seq_len(n)), 1, paste, collapse = '-') ## TODO: verify order and the similar approaches above
    ait <- apply(expand.grid(seq_len(n), seq_len(nt)), 1, paste, collapse = '-')

    t <- as.numeric(factor(t))
    i <- as.numeric(i)
    j <- as.numeric(j)

    res <- sapply(1:length(value), function(index) {
        wij <- which(paste(i[index], j[index], sep = '-') == aij)
        wit <- which(paste(i[index], t[index], sep = '-') == ait)
        wjt <- which(paste(j[index], t[index], sep = '-') == ait)
        Tij <- length(which(i == i[index] & j == j[index]))
        value[index] -
            mean(value[which(i == i[index] & j == j[index])]) -
                (Ab[, wij] %*% fib) / Tij - fib[wit] - o[wjt] +
                    (Matrix::t(A[, wij]) %*% o) / Tij + k[wit, ] -
                        (Matrix::t(Ab[, wij]) %*% k) / Tij
    })

    res

}
