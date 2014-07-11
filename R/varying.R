#' Varying model based on formula (6)
#' @param mx a \code{matrix} holding 4 columns for \code{i}, \code{j}, \code{t} and \code{value}
#' @param my a \code{matrix} holding 4 columns for \code{i}, \code{j}, \code{t} and \code{value}
#' @export
#' @return numeric value of OLS beta
varyingModel1 <- function(mx, my) {

    ## TODO: generalize to work if K>1

    l <- levels(mx[, 1])
    N <- length(l)
    T <- length(unique(mx[, 3]))

    ## sort data by i, j, t
    mx <- mx[order(mx[, 1], mx[, 2], mx[, 3]), ]

    Aij <- lapply(structure(l, .Names = l), function(i) {
        lapply(structure(l, .Names = l), function(j) {
            sum(sapply(1:T, function(t) {
                w <- which(mx[, 1] == i & mx[, 2] == j & mx[, 3] == t)
                mx[w, 4] * t(mx[w, 4])
            }))
        })
    }) ## N x N
    Cxx <- Reduce('+', lapply(unlist(Aij, recursive = FALSE), solve))

    ## MX <- N^2 %*% t(mx[, 4]) %*% solve(Aij) %*% solve(Cxx)
    MX <- as.vector(sapply(l, function(i) {
        sapply(l, function(j) {
            sapply(1:T, function(t) {
                w <- which(mx[, 1] == i & mx[, 2] == j & mx[, 3] == t)
                N^2 %*% t(mx[w, 4]) %*% solve(Aij[[i]][[j]]) %*% solve(Cxx)
            })})})) ## TODO check order

}
