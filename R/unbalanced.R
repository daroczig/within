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
#' @aliases WithinTransformation1c WithinTransformation2c WithinTransformation3c WithinTransformation4c WithinTransformation5c WithinTransformation6c
WithinTransformation1c <- function(i = m[, 1], j = m[, 2], t = m[, 3], value = m[, 4:ncol(m)], m)  { # formula (22)

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

    l <- sapply(Vs, function(x) x@Dim[1])
    names(l) <- 1:length(l)
    js <- as.numeric(rep(names(l), times = l))

    D2 <- sparseMatrix(1:length(js), js, x = 1)

    dN2 <- t(D1) %*% D1

    ## inverse of dN2 with reduced number of cells
    ix <- union(which(rowSums(dN2) > 0), which(colSums(dN2) > 0))
    dN2i <- dN2
    dN2i[ix, ix] <- solve(dN2[ix, ix])

    Aa  <- t(D2) %*% D1
    Da <- D2 - (D1 %*% dN2i %*% t(Aa))


    dTx <- t(D2) %*% D2
    It  <- sparseMatrix(1:dim(D1)[1], 1:dim(D1)[1], x = 1)
    Qa <- t(D2) %*% Da

    ## geninv of Qa
    Qai <- ginv(as.matrix(Qa))
    Qai <- Matrix(Qai)

    ## full value (added NA)
    fulli <- rep(NA, )

    ## this should work on smaller matrices (formula 22)
    ## Pa <- (It - (D1 %*% dN2i %*% t(D1))) - (Da %*% Qai %*% t(Da))


    fia <- Qai %*% t(Da) %*% value

    t <- as.numeric(factor(t))
    i <- as.numeric(i)
    j <- as.numeric(j)

    res <- sapply(1:length(value), function(index) {
        value[index] - mean(value[which(i == i[index] & j == j[index])]) - fia[as.numeric(t[index])]
    })

}
