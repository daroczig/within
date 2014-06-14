#' Within transformation with missing data
#'
#' Either \code{i}, \code{j}, \code{t} and \code{value} parameters or \code{m} data frame with all these four vectors are required.
#' @param i country 1 (vector)
#' @param j country 2 (vector)
#' @param t year (vector)
#' @param value vector to be transformed
#' @param m a \code{matrix} holding 4 columns for \code{i}, \code{j}, \code{t} and \code{value}
#' @return vector
#' @references László Balázsi, László Mátyás and Tom Wansbeek (2014): The Estimation of Multi-dimensional Fixed Effects Panel Data Models. Formula (19), (20), and (21).
#' @export
#' @aliases WithinTransformation1b WithinTransformation2b WithinTransformation3b WithinTransformation4b WithinTransformation5b WithinTransformation6b
WithinTransformation1b <- function(i = m[, 1], j = m[, 2], t = m[, 3], value = m[, 4:ncol(m)], m) # formula (19)
    value - (N()-1) / (N() * (N()-2)*T()) * (`yi++`()+`y+j+`()) - 1 / (N() * (N()-2)*T()) * (`yj++`() + `y+i+`()) - 1 / (N() * (N()-1)) * (`y++t`()) + 2 / (N() * (N()-2)*T()) * (`y+++`())
#' @export
WithinTransformation2b <- WithinTransformation2
#' @export
WithinTransformation3b <- WithinTransformation3
#' @export
WithinTransformation4b <- WithinTransformation4
#' @export
WithinTransformation5b <- function(i = m[, 1], j = m[, 2], t = m[, 3], value = m[, 4:ncol(m)], m) # formula (20)
    value - (N()-1) / (N() * (N()-2)) * (`yi+t`()+`y+jt`()) - 1 / (N() * (N()-2)) * (`y+it`() + `yj+t`()) + 1 / ((N()-1) * (N()-2)) * (`y++t`())
#' @export
WithinTransformation6b <- function(i = m[, 1], j = m[, 2], t = m[, 3], value = m[, 4:ncol(m)], m) # formula (21)
    value - (1/(N()-1))*`y+jt`() - (1/(N()-1))*`yi+t`() - 1/T()*`yij+`() + (1/((N()-1)^2))*`y++t`() + (1/((N()-1)*T()))*`y+j+`() + (1/((N()-1)*T()^2))*`y+++`() + (1/((N()-1)*T()))*`yji+`() - (1/(N()-1))*`yjit`()

## helper functions
`T` <- function() eval.parent(expression(length(table(t))))
`N` <- function() eval.parent(expression(length(levels(i))))
`yi++` <- function() eval.parent(expression(ave(value, i, FUN = sum, na.rm = TRUE)))
`y+j+` <- function() eval.parent(expression(ave(value, j, FUN = sum, na.rm = TRUE)))
`yj++` <- function() eval.parent(expression({
    m$id <- 1:nrow(m)
    a <- aggregate(value ~ i, FUN = sum)
    names(a) <- c('j', 'res')
    ma <- merge(m, a)
    ma$res[order(ma$id)]
}))
`y+i+` <- function() eval.parent(expression({
    m$id <- 1:nrow(m)
    a <- aggregate(value ~ j, FUN = sum)
    names(a) <- c('i', 'res')
    ma <- merge(m, a)
    ma$res[order(ma$id)]
}))
`y++t` <- function() eval.parent(expression(ave(value, t, FUN = sum, na.rm = TRUE)))
`y+++` <- function() eval.parent(expression(sum(value, na.rm = TRUE)))
`yi+t` <- function() eval.parent(expression(ave(value, i, t, FUN = sum, na.rm = TRUE)))
`y+jt` <- function() eval.parent(expression(ave(value, j, t, FUN = sum, na.rm = TRUE)))
`y+it` <- function() eval.parent(expression({
    m$id <- 1:nrow(m)
    a <- aggregate(value ~ j + t, FUN = sum)
    names(a) <- c('i', 't', 'res')
    ma <- merge(m, a)
    ma$res[order(ma$id)]
}))
`yj+t` <- function() eval.parent(expression({
    m$id <- 1:nrow(m)
    a <- aggregate(value ~ i + t, FUN = sum)
    names(a) <- c('j', 't', 'res')
    ma <- merge(m, a)
    ma$res[order(ma$id)]
}))
`yjit` <- function() eval.parent(expression(value))
`yij+` <- function() eval.parent(expression(ave(value, i, j, FUN = sum, na.rm = TRUE)))
`yji+` <- function() eval.parent(expression({
    m$id <- 1:nrow(m)
    a <- aggregate(value ~ i + j, FUN = sum)
    names(a) <- c('j', 'i', 'res')
    ma <- merge(m, a)
    ma$res[order(ma$id)]
}))
