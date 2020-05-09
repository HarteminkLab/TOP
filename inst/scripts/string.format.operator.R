
# define string concatenation operator
`%^%` <- function(fmt, ...) {
    do.call(sprintf, as.list(c(fmt , ...)))
}

concat <- function(..., sep='') {
  paste(..., sep=sep, collapse=sep)
}
