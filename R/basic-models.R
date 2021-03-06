#' OLS on transformed vectors
#' @param mx a \code{matrix} holding at least 4 columns for \code{i}, \code{j}, \code{t} and any number of \code{value}(s)
#' @param my a \code{matrix} holding 4 columns for \code{i}, \code{j}, \code{t} and \code{value}
#' @param transformation tranforming function, e.g. \code{WithinTransformation2}
#' @return numeric (beta)
#' @export
#' @examples \dontrun{
#' mx <- cbind(m$country1, m$country2, m$year, m$gattwto1)
#' my <- cbind(m$country1, m$country2, m$year, m$lnrgdp1)
#' OLSonTransformation(mx, my, WithinTransformation2)
#' OLSonTransformation(my, mx, WithinTransformation2)
#' }
#' @importFrom digest digest
OLSonTransformation <- function(mx, my, transformation, checkpointing = FALSE) {

    mxn <- tail(names(mx), -3)

    ## transform
    my <- transformation(m = my)
    pd <- digest(list(my, transformation), algo = 'sha1')
    if (isTRUE(checkpointing))
        saveRDS(my, paste0('my-', pd))
    if (require(parallel) & ncol(mx) > 4) {
        mx <- simplify2array(mclapply(4:ncol(mx), function(i) {
            r <- transformation(m = mx[, c(1:3, i)])
            pd <- digest(list(mx[, c(1:3, i)], transformation), algo = 'sha1')
            if (isTRUE(checkpointing))
                saveRDS(r, paste0('mx-', names(mx)[i], '-', pd))
            r
        }))
    } else {
        mx <- sapply(4:ncol(mx), function(i) transformation(m = mx[, c(1:3, i)]))
    }

    ## helper vars
    nval <- ncol(mx)
    MX <- mx

    ## remove transformed vars with zero variance
    na <- which(apply(mx, 2, var) == 0)
    if (length(na) > 0) {
        mx <- mx[, -na]
    }

    ## beta
    b <- solve(t(mx) %*% as.matrix(mx)) %*% t(mx) %*% as.matrix(my)

    ## standard error
    se <- as.numeric(t(my - mx %*% b) %*% (my - mx %*% b)) * solve(t(mx) %*% mx) / (nrow(mx) - ncol(mx))
    se <- sqrt(diag(se))

    ## R squared
    r2 <- t(mx %*% b) %*% (mx %*% b) / t(my) %*% my

    ## add NA beta coefficient for variables with zero variance
    reAddNA <- function(v) {
        t <- matrix(NA, nval, 1)
        t[setdiff(1:nval, na), 1] <- v
        t
    }
    if (length(na) > 0) {
        b  <- reAddNA(b)
        se <- reAddNA(se)
    }

    ## return
    data.frame(beta = b, se = se, r2 = r2, row.names = mxn)

}


#' Within transformation
#'
#' Either \code{i}, \code{j}, \code{t} and \code{value} parameters or \code{m} data frame with all these four vectors are required.
#' @param i country 1 (vector)
#' @param j country 2 (vector)
#' @param t year (vector)
#' @param value vector to be transformed
#' @param m a \code{matrix} holding 4 columns for \code{i}, \code{j}, \code{t} and \code{value}
#' @return vector
#' @references László Balázsi, László Mátyás and Tom Wansbeek (2014): The Estimation of Multi-dimensional Fixed Effects Panel Data Models. Formula (3), (6), (9), (11) and (13).
#' @examples \dontrun{
#' str(WithinTransformation2(m=cbind(m$country1, m$country2, m$year, m$lnrgdp1)))
#' str(WithinTransformation2(m$country1, m$country2, m$year, m$lnrgdp1))
#' }
#' @export
#' @aliases WithinTransformation1 WithinTransformation2 WithinTransformation3 WithinTransformation4 WithinTransformation5 WithinTransformation6
WithinTransformation1 <- function(i = m[, 1], j = m[, 2], t = m[, 3], value = m[, 4:ncol(m)], m) # model (1) formula (3)
    value - yi..() - y.j.() - y..t() + 2*y...()
#' @export
WithinTransformation2 <- function(i = m[, 1], j = m[, 2], t = m[, 3], value = m[, 4:ncol(m)], m) # model (5) formula (6)
    value - yij.()
#' @export
WithinTransformation3 <- function(i = m[, 1], j = m[, 2], t = m[, 3], value = m[, 4:ncol(m)], m) # model (7) formula ()
    value - yij.() - y..t() + y...()
#' @export
WithinTransformation4 <- function(i = m[, 1], j = m[, 2], t = m[, 3], value = m[, 4:ncol(m)], m) # model (8) formula (9)
    value - yi.t()
#' @export
WithinTransformation5 <- function(i = m[, 1], j = m[, 2], t = m[, 3], value = m[, 4:ncol(m)], m) # model (10) formula (11)
    value - yi.t() - y.jt() + y..t()
#' @export
WithinTransformation6 <- function(i = m[, 1], j = m[, 2], t = m[, 3], value = m[, 4:ncol(m)], m) # model (12) formula (13)
    value - yij.() - yi.t() - y.jt() + yi..() + y.j.() + y..t() - y...()


## helper functions
`yi..` <- function() eval.parent(expression(ave(value, i, FUN = mean, na.rm = TRUE)))
`y.j.` <- function() eval.parent(expression(ave(value, j, FUN = mean, na.rm = TRUE)))
`y..t` <- function() eval.parent(expression(ave(value, t, FUN = mean, na.rm = TRUE)))
`yij.` <- function() eval.parent(expression(ave(value, i, j, FUN = mean, na.rm = TRUE)))
`yi.t` <- function() eval.parent(expression(ave(value, i, t, FUN = mean, na.rm = TRUE)))
`y.jt` <- function() eval.parent(expression(ave(value, j, t, FUN = mean, na.rm = TRUE)))
`y...` <- function() eval.parent(expression(mean(value, na.rm = TRUE)))
