#' Varying model based on formula (6)
#' @param mx a \code{matrix} holding 4 columns for \code{i}, \code{j}, \code{t} and \code{value}
#' @param my a \code{matrix} holding 4 columns for \code{i}, \code{j}, \code{t} and \code{value}
#' @export
#' @return numeric value of OLS beta
#' @references László Balázsi, László Mátyás (2014): The Estimation of Varying Coefficients Multi-dimensional Panel Data Models. Formula (6).
varyingModel1 <- function(mx, my) {

    if (nrow(mx) != nrow(my))
        stop('incompatible matrices provided')
    if (!identical(mx[, 1:3], my[, 1:3]))
        stop('incompatible matrices provided')

    l <- levels(mx[, 1])
    N <- length(l)
    T <- length(unique(mx[, 3]))
    K <- ncol(mx) - 3

    ## sort data by i, j, t
    mx <- mx[order(mx[, 1], mx[, 2], mx[, 3]), ]

    Aij <- lapply(structure(l, .Names = l), function(i) {
        lapply(structure(l, .Names = l), function(j) {
            sum(sapply(1:T, function(t) {
                w <- which(mx[, 1] == i & mx[, 2] == j & mx[, 3] == t)
                mx[w, 3+(1:K)] * t(mx[w, 3+(1:K)])
            }))
        })
    }) ## N x N
    Cxx <- Reduce('+', lapply(unlist(Aij, recursive = FALSE), solve))

    MX <- as.vector(sapply(l, function(i) {
        sapply(l, function(j) {
            sapply(1:T, function(t) {
                w <- which(mx[, 1] == i & mx[, 2] == j & mx[, 3] == t)
                N^2 %*% t(mx[w, 3+(1:K)]) %*% solve(Aij[[i]][[j]]) %*% solve(Cxx)
            })})})) ## TODO check order

    Bij <- lapply(structure(l, .Names = l), function(i) {
        lapply(structure(l, .Names = l), function(j) {
            sum(sapply(1:T, function(t) {
                w <- which(mx[, 1] == i & mx[, 2] == j & mx[, 3] == t)
                mx[w, 3+(1:K)] * my[w, 4]
            }))
        })
    })
    Cxy <- sum(mapply(function(x, y) solve(x) %*% y, unlist(Aij, recursive = FALSE), unlist(Bij, recursive = FALSE)))

    MY <- as.vector(sapply(l, function(i) {
        sapply(l, function(j) {
            sapply(1:T, function(t) {
                w <- which(mx[, 1] == i & mx[, 2] == j & mx[, 3] == t)
                my[w, 4] - t(mx[w, 3+(1:K)]) %*% ((solve(Aij[[i]][[j]]) %*% Bij[[i]][[j]]) - (solve(Aij[[i]][[j]]) %*% solve(Cxx) %*% Cxy))
            })})})) ## TODO check order

    ## OLS
    solve(t(MX) %*% as.matrix(MX)) %*% t(MX) %*% as.matrix(MY)

}


#' Varying model based on formula (12)
#' @param mx a \code{matrix} holding 4 columns for \code{i}, \code{j}, \code{t} and \code{value}
#' @param my a \code{matrix} holding 4 columns for \code{i}, \code{j}, \code{t} and \code{value}
#' @export
#' @return OLS beta
#' @references László Balázsi, László Mátyás (2014): The Estimation of Varying Coefficients Multi-dimensional Panel Data Models. Formula (12).
varyingModel3 <- function(mx, my, mz) {

    if (nrow(mx) != nrow(my) | nrow(mx) != nrow(mz) | nrow(my) != nrow(mz))
        stop('incompatible matrices provided')
    if (!identical(mx[, 1:3], my[, 1:3]) | !identical(mx[, 1:3], mz[, 1:3]) | !identical(my[, 1:3], mz[, 1:3]))
        stop('incompatible matrices provided')

    l <- levels(mx[, 1])
    N <- length(l)
    T <- length(unique(mx[, 3]))
    K <- ncol(mx) - 3

    ## sort data by i, j, t
    mx <- mx[order(mx[, 1], mx[, 2], mx[, 3]), ]

    Aij <- lapply(structure(l, .Names = l), function(i) {
        lapply(structure(l, .Names = l), function(j) {
            sum(sapply(1:T, function(t) {
                w <- which(mx[, 1] == i & mx[, 2] == j & mx[, 3] == t)
                mx[w, 3+(1:K)] * t(mx[w, 3+(1:K)])
            }))
        })
    }) ## N x N
    Cxx <- Reduce('+', lapply(unlist(Aij, recursive = FALSE), solve))

    MX <- as.vector(sapply(l, function(i) {
        sapply(l, function(j) {
            sapply(1:T, function(t) {
                w <- which(mx[, 1] == i & mx[, 2] == j & mx[, 3] == t)
                N^2 %*% t(mx[w, 3+(1:K)]) %*% solve(Aij[[i]][[j]]) %*% solve(Cxx)
            })})})) ## TODO check order

    Bij <- lapply(structure(l, .Names = l), function(i) {
        lapply(structure(l, .Names = l), function(j) {
            sum(sapply(1:T, function(t) {
                w <- which(mx[, 1] == i & mx[, 2] == j & mx[, 3] == t)
                mx[w, 3+(1:K)] * my[w, 4]
            }))
        })
    })
    Cxy <- sum(mapply(function(x, y) solve(x) %*% y, unlist(Aij, recursive = FALSE), unlist(Bij, recursive = FALSE)))

    MY <- as.vector(sapply(l, function(i) {
        sapply(l, function(j) {
            sapply(1:T, function(t) {
                w <- which(mx[, 1] == i & mx[, 2] == j & mx[, 3] == t)
                my[w, 4] - t(mx[w, 3+(1:K)]) %*% ((solve(Aij[[i]][[j]]) %*% Bij[[i]][[j]]) - (solve(Aij[[i]][[j]]) %*% solve(Cxx) %*% Cxy))
            })})})) ## TODO check order

    Dij <- lapply(structure(l, .Names = l), function(i) {
        lapply(structure(l, .Names = l), function(j) {
            sum(sapply(1:T, function(t) {
                w <- which(mx[, 1] == i & mx[, 2] == j & mx[, 3] == t)
                mx[w, 3+(1:K)] * t(mz[w, 4])
            }))
        })
    })
    Cxz <- sum(mapply(function(x, y) solve(x) %*% y, unlist(Aij, recursive = FALSE), unlist(Dij, recursive = FALSE)))

    MZ <- as.vector(sapply(l, function(i) {
        sapply(l, function(j) {
            sapply(1:T, function(t) {
                w <- which(mx[, 1] == i & mx[, 2] == j & mx[, 3] == t)
                t(mz[w, 4]) - t(mx[w, 3+(1:K)]) %*% ((solve(Aij[[i]][[j]]) %*% Dij[[i]][[j]]) - (solve(Aij[[i]][[j]]) %*% solve(Cxx) %*% Cxz))
            })})})) ## TODO check order

    ## OLS

}
