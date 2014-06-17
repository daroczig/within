#' Within transformation for dynamic autoregressive models
#'
#' Either \code{i}, \code{j}, \code{t} and \code{value} parameters or \code{m} data frame with all these four vectors are required.
#' @param i country 1 (vector)
#' @param j country 2 (vector)
#' @param t year (vector)
#' @param value vector to be transformed
#' @param m a \code{matrix} holding 4 columns for \code{i}, \code{j}, \code{t} and \code{value}
#' @return vector
#' @references László Balázsi, László Mátyás and Tom Wansbeek (2014): The Estimation of Multi-dimensional Fixed Effects Panel Data Models. pp. 22--24.
#' @export
#' @examples {
#' m <- data.frame(
#'     i     = rep(letters[1:3], each = 3*5),
#'     j     = rep(letters[1:3], times = 3*5),
#'     t     = rep(1:5, each = 3),
#'     value = runif(3*3*5)*10)
#' }
#' @aliases WithinTransformation1d WithinTransformation2d WithinTransformation3d WithinTransformation4d WithinTransformation5d WithinTransformation6d
WithinTransformation1d <- WithinTransformation1
#' @export
WithinTransformation2d <- function(i = m[, 1], j = m[, 2], t = m[, 3], value = m[, 4:ncol(m)], m) {

    ## global vars
    n  <- length(unique(i))
    nt <- length(table(t))

    ## check if provided matrix is valid
    if (length(unique(i)) != length(unique(j)))
        stop('i and j categories differ.')

    ## check if we have all data for all t
    missing_table <- table(t, paste(i, j, sep = '_'))
    missing_pairs <- dimnames(missing_table)[[2]][unique(which(missing_table == 0, arr.ind = TRUE)[, 2])]
    w <- which(paste(i, j, sep = '_') %in% missing_pairs)
    if (length(w) > 1) {
        warning(paste('Some missing "i"-"j" pairs identified a few "t", thus removing', length(w), 'records'))
        for (v in c('i', 'j', 't', 'value'))
            assign(v, get(v)[-w])
    }

    ## reorder by t, i and j
    o <- order(t, i, j)
    for (v in c('i', 'j', 't', 'value'))
        assign(v, get(v)[o])

    ## deltaY
    M <- data.frame(
        i = i,
        j = j,
        t = t,
        value = value
        )
    M <- M[order(M$i, M$j), ]
    M$d1 <- unlist(by(M$value, list(M$i, M$j) , function(x) c(NA, diff(x))), use.names = FALSE)
    M$d2 <- unlist(by(M$d1, list(M$i, M$j) , function(x) c(NA, diff(x))), use.names = FALSE)

    dY1 <- head(as.numeric(na.omit(M$d1)), -1*n^2)
    dY2 <- as.numeric(na.omit(M$d2))
    rm(M)

    ## z and Z
    z <- Matrix(0, nrow = nt-2, ncol = (nt-1)*(nt-2)/2)
    a <- apply(expand.grid(unique(i), unique(i)), 1, paste, collapse = '-')
    Zs <- lapply(a, function(w) {
        wi <- strsplit(w, '-')[[1]][1]
        wj <- strsplit(w, '-')[[1]][2]
        w  <- tail(which(as.character(i) == wi & as.character(j) == wj), -2)

        if (length(w) != (nt-2))
            stop(paste0('Missing data for "i=', wi, '" and "j=', wj, '" pairs.'))

        z <- Matrix(0, nrow = nt-2, ncol = (nt-1)*(nt-2)/2)
        for (ti in seq_len(nt-2)) {
            to   <- tail(cumsum(1:ti), 1)
            from <- to - ti + 1
            z[ti, from:to] <- head(value[w], n = ti)
        }

        return(z)
    })
    Z <- do.call(rBind, Zs)

    ## Σ matrix
    S <- matrix(0, nrow = nt-2, ncol = nt-2)
    diag(S) <- 2
    diag(S[-nrow(S),-1]) <- diag(S[-1,-ncol(S)]) <- -1

    ## Ω matrix
    O <- kronecker(sparseMatrix(1:n^2, 1:n^2, x = 1), S)

    solve(t(dY2) %*% Z %*% solve(t(Z) %*% O %*% Z) %*% t(Z) %*% dY2) %*%
        t(dY2) %*% Z %*% solve(t(Z) %*% O %*% Z) %*% t(Z) %*% dY1

}
#' @export
WithinTransformation3d <- function(i = m[, 1], j = m[, 2], t = m[, 3], value = m[, 4:ncol(m)], m) {



}
#' @export
WithinTransformation4d <- WithinTransformation4
#' @export
WithinTransformation5d <- WithinTransformation5
#' @export
WithinTransformation6d <- function(i = m[, 1], j = m[, 2], t = m[, 3], value = m[, 4:ncol(m)], m) {



}
